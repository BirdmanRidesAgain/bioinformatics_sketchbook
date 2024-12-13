#!/bin/bash
# Keiler Collier
# 30 July 2024
# This is a shell script designed to run mapping on paired-end fastq files.
# It assumes that your conda environment can find all the tools necessary, and that all fasta inputs are in the directory you hand in
# it also 

# If no arguments are given, give usage and exit
USAGE="Usage: bash map.sh <sample> <refseq.fq> <read1> <read2>"
OUTPUTS="outputs .bam and .in directory <sample>"
if [ -z $1 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    exit 1
fi

####################################################################################
### SETTING VARIABLES
SAMPLE=$1
REFSEQ=$2
READ1=$3
READ2=$4

THREADS=80


####################################################################################
### SOURCING CONDA
eval "$(conda shell.bash hook)"

####################################################################################
### TRIM WITH FASTP

if [ ! -f ${SAMPLE}/fastq-files/fastp.json ]
then 
    conda activate fastp

    # Known bug - if you activate fastp on a file that already has a fastp'd name, it will erase the file.
    fastp -i $READ1 -I $READ2 -o ${SAMPLE}/${SAMPLE}_R1_fastp.fq.gz -O ${SAMPLE}/${SAMPLE}_R2_fastp.fq.gz
    READ1=${SAMPLE}/fastq-files/${SAMPLE}_R1_fastp.fq.gz
    READ2=${SAMPLE}/fastq-files/${SAMPLE}_R2_fastp.fq.gz
    mv fastp* ${SAMPLE}/fastq-files

    conda deactivate
fi

########################################################################
### BAM PROCESSING WITH SAMTOOLS ONLY (NO PICARD, NO GATK):
FINAL_BAM=${SAMPLE}/${SAMPLE}_final.bam

if [ ! -f ${FINAL_BAM} ]
then
    conda activate genome_assembly

    echo "Mapping ${SAMPLE} with bwa-mem"
    bwa mem -t $THREADS $REFSEQ $READ1 $READ2 | samtools view -@ $THREADS -bS > ${SAMPLE}/${SAMPLE}.bam

    echo "Cleaning ${SAMPLE} bam with samtools"
    samtools collate -@ $THREADS -o ${SAMPLE}/${SAMPLE}_namecollate.bam ${SAMPLE}/${SAMPLE}.bam
    samtools fixmate -@ $THREADS -m ${SAMPLE}/${SAMPLE}_namecollate.bam ${SAMPLE}/${SAMPLE}_fixmate.bam
    samtools sort -@ $THREADS -o ${SAMPLE}/${SAMPLE}_positionsort.bam ${SAMPLE}/${SAMPLE}_fixmate.bam
    samtools markdup -r -@ $THREADS ${SAMPLE}/${SAMPLE}_positionsort.bam ${SAMPLE}/${SAMPLE}_markdup.bam
    samtools addreplacerg -@ $THREADS -r ID:${SAMPLE} -r LB:${SAMPLE} -r SM:${SAMPLE} -o ${FINAL_BAM} ${SAMPLE}/${SAMPLE}_markdup.bam

    # Do flagstat for the hell of it
    samtools flagstat ${FINAL_BAM} > ${SAMPLE}/${SAMPLE}_final_samtools_flagstat.txt

    # Removing intermediate bamfiles:
    rm ${SAMPLE}/${SAMPLE}_{namecollate,fixmate,positionsort,markdup}.bam 
    rm ${SAMPLE}/${SAMPLE}.bam

    conda deactivate
fi

#### GETTING METRICS ON THE FINAL BAM:
BAMQC=${SAMPLE}/${SAMPLE}_bamqc/genome_results.txt
if [ ! -f ${BAMQC} ]
then
    conda activate qualimap

    echo "Evaluating ${FINAL_BAM} with bamqc"
    qualimap bamqc -bam ${FINAL_BAM} -outdir ${SAMPLE}/${SAMPLE}_bamqc -outformat PDF:HTML --java-mem-size=10G

    conda deactivate
fi
echo "Mapping complete for ${SAMPLE}"
