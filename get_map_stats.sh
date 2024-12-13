#!/bin/bash
# Keiler Collier
# 13 Dec 2024
# get_map_stats.sh
# This is a shell script designed to compile summary statistic results from many bamfiles.
# It depends on the outputs from fastp and bamqc, which we get from the 'map.sh' script.
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


USAGE="Usage: get_map_stats.sh <output_prefix> <refseq.fq> <samplelist>."
OUTPUTS="Outputs oneline summary tsv."
if [ -z $1 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

####################################################################################
### SETTING VARIABLES
OUTPUT_PREFIX=$1
REFSEQ=$2
SAMPLELIST=$3

if [ -z $REFSEQ]
then
    echo "No refseq supplied."
    exit 1
fi

# Get a list of all the samples, if it was not provided as an argument
if [ -z $SAMPLELIST ]
then
    echo -e "No samplelist provided. Generating.\n"
    find . -maxdepth 1 -type d  -printf "%f\n" | grep -v "bamqc" | sort | tail -n +2 > ${OUTPUT_PREFIX}_samplelist.txt
    SAMPLELIST=${OUTPUT_PREFIX}_samplelist.txt
fi

########################################################################
HEADER="bamfile\trefseq\tdepth\tdepth_stdev\tnum_reads\tpercent_mapped\tpercent_duplicate\tavg_mapqual"
# Output header for our summary file
echo -e $HEADER > ${OUTPUT_PREFIX}_mapping_stats.tsv;

while read SAMPLE; do
    #### GETTING METRICS ON THE FINAL BAM:
    echo -e $HEADER > ${SAMPLE}/${SAMPLE}_mapping_stats.tsv;

    # Check if the statfiles exist
    BAMQC=${SAMPLE}/${SAMPLE}_bamqc/genome_results.txt;
    FASTP=${SAMPLE}/fastq-files/fastp.json

    # figure out which one (or both) were missing
    if [ ! -f $BAMQC ] || [ ! -f $FASTP ]
    then
        echo -en "${SAMPLE} qualstats_error\t" >> ${SAMPLE}/${SAMPLE}_mapping_stats.tsv;
        if [ -f $BAMQC ]
        then
            echo -e "Fastp missing" >> ${SAMPLE}/${SAMPLE}_mapping_stats.tsv;
        elif [ -f $FASTP ]
        then
            echo -e "Bamqc missing" >> ${SAMPLE}/${SAMPLE}_mapping_stats.tsv;
        else
            echo -e "Both bamqc and fastp missing" >> ${SAMPLE}/${SAMPLE}_mapping_stats.tsv;
        fi
    else
        #### FORMAT METRICS INTO TSV
        REFSEQ=$(basename $REFSEQ)
        DEPTH=$(cat ${BAMQC} | grep "mean coverageData" | awk '{print $NF}' | sed s'/X//');
        DEPTH_STDEV=$(cat ${BAMQC} | grep "std coverageData" | awk '{print $NF}' | sed s'/X//');
        NUM_RDS=$(cat $BAMQC | grep "number of reads" | sed s'/,//'g | awk '{print $NF}');
        PERCENT_PRIMARY_MAPPED=$(cat ${BAMQC} | grep "number of mapped reads" | awk '{print $NF}' | sed s'/(//' | sed s'/)//' | sed s'/%//');
        PRIMARY_DUPLICATES=$(grep "rate" ${FASTP} | grep -v "q" | awk '{print $NF}');

        AVG_MAPQUAL=$(cat ${BAMQC} | grep "mean mapping quality" | awk '{print $NF}');
        echo -e "${SAMPLE}_final.bam\t${REFSEQ}\t${DEPTH}\t${DEPTH_STDEV}\t${NUM_RDS}\t${PERCENT_PRIMARY_MAPPED}\t${PRIMARY_DUPLICATES}\t${AVG_MAPQUAL}" >> ${SAMPLE}/${SAMPLE}_mapping_stats.tsv;
    fi
done < $SAMPLELIST

# Now cat all the information into the master file
for i in */*_mapping_stats.tsv; do
    tail $i -n +2 >> ${OUTPUT_PREFIX}_mapping_stats.tsv;
done
