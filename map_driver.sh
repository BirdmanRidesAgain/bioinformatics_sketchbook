#!/bin/bash
# Keiler Collier
# 13 Dec 2024
# map_driver.sh
# This is a shell script designed to 'drive' map.sh and get_map_stats.sh to partially automate our mapping.
# You activate it in a directory that follows the appropriate pattern, and it tries to find things to map.
#
# You are expected to supply a samplelist - this can be obtained with:
#   "find . -type d -maxdepth 1 -printf "%f\n" | grep -v "bamqc" | sort | tail -n +2 > samplelist"
# 
# The samplelist should consist of one individual per line:
# ind1
# ind2
# ind3
#
# The following directory structure is expected:
#   ${SAMPLE}
#       ${SAMPLE}_final.bam
#       ${SAMPLE}_bamqc
#           genome_results.txt
#       fastq-files
#           fastp.json
#           ${SAMPLE}*.fq.gz
#

USAGE="Usage: map_driver.sh <output_prefix> <refseq.fq> <samplelist>."
OUTPUTS="Outputs oneline summary tsv."
if [ -z $1 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

####################################################################################
### SETTING VARIABLES
OUTPUT_PREFIX=$1;
REFSEQ=$2;
SAMPLELIST=$3;

####################################################################################
### MAIN FUNCTION
# Get a list of all the samples, if it was not provided as an argument
if [ -z $SAMPLELIST ]
then
    echo -e "No samplelist provided. Generating.\n";
    find . -maxdepth 1 -type d  -printf "%f\n" | grep -v "bamqc" | sort | tail -n +2 > ${OUTPUT_PREFIX}_samplelist.txt;
    SAMPLELIST=${OUTPUT_PREFIX}_samplelist.txt;
fi;

####################################################################################
### ATTEMPT TO MAP
# Get a list of all the fq files per sample (sometimes weird stuff happens in data processing and we want only paired reads)
while read SAMPLE; do 
    echo -en "${SAMPLE} "
    echo $(find ./${SAMPLE}/fastq-files/* -name *fq.gz);
done < $SAMPLELIST > ${OUTPUT_PREFIX}_fqlist.txt
FQLIST=${OUTPUT_PREFIX}_fqlist.txt

# We will sanity-check the number (and matching) of FQ files before attempting to map
while read FQs; do
    NUM_FQ_FILES=$(echo $FQs | awk '{print NF -1}')
    SAMPLE=$(echo $FQs | awk '{print $1}')
    if [ $NUM_FQ_FILES -eq 2 ] # we assume that the two fq files are properly paired. This will fail if they aren't.
    then
        if [ ! -f ${SAMPLE}/${SAMPLE}_final.bam ]
        then
            echo "Mapping ${SAMPLE}"
            READ1=$($FQs | awk '{print $2}')
            READ2=$($FQs | awk '{print $3}')
            map.sh ${SAMPLE} ${REFSEQ} ${READ1} ${READ2}
        fi
    else 
        echo "${SAMPLE} does not have two .fq files. Check it for quality."
    fi
done < $FQLIST

####################################################################################
### ATTEMPT TO GET VARSTATS
get_map_stats.sh ${OUTPUT_PREFIX} ${REFSEQ} ${SAMPLELIST}
echo "end script."