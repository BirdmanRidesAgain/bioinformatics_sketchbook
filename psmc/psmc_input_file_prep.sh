#!/bin/bash
# usage: psmc_input_file_prep.sh <assembly.fa> <reads.fq.gz> <project_name>
# outputs:

########################################################################
### SETTING VARIABLES
# Positial variables
ASSEMBLY=$1
READS=$2
PROJECT_NAME=$3

# Hardcoded variables
THREADS=40
LOWER_BOUND=10
UPPER_BOUND=100


minimap2 -t $THREADS -x map-ont $ASSEMBLY $READS | samtools sort | samtools view -bS > ${PROJECT_NAME}_ONT_remapped.bam
MAPFILE=${PROJECT_NAME}_ONT_remapped.bam
bcftools mpileup -C50 -uf $ASSEMBLY $MAPFILE | bcftools view -c - | vcfutils.pl vcf2fq -d $LOWER_BOUND -D $UPPER_BOUND | gzip > ${PROJECT_NAME}.fq.gz

### Command used for experimentation purposes
#bcftools mpileup -Ou -f $ASSEMBLY -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR $MAPFILE | bcftools call -mv -Oz -f GQ,GP -o Burhinus_ONT_remapped_sort.vcf.gz
