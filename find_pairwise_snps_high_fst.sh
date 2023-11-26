#!/bin/bash
# usage: find_pairwise_snps_high_fst.sh <input.vcf> <population1.txt> <population2.txt>
# outputs: <${COMPARISON}_fixed-SNPs.txt> <${COMPARISON}_high_effect-SNPs.txt>

# Script currently finds Fst values with a pairwise fixation greater than 70%. This is hardcoded.
####################################################################################
### SOURCING VARIABLES AND CONDA
source  ~/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

########################################################################
### SETTING VARIABLES
VCF_FILE=$1
POP_FILE1=$2
POP_FILE2=$3
MIN_FST_VALUE=7

### PROCESSING VARIABLES
POP1=`basename $POP_FILE1 | sed 's/\.txt//g'`
POP2=`basename $POP_FILE2 | sed 's/\.txt//g'`
COMPARISON=${POP1}_${POP2}

########################################################################
### ECHOING PARAMETERS
echo "Parameters as interpreted:"
echo ""
echo "Input vcf: $VCF_FILE"
echo "Population 1 from: $POP_FILE1"
echo "Population 2 from: $POP_FILE2"
echo "Output name: $COMPARISON"
echo ""

########################################################################
### MAKE DIRECTORIES
OUT_DIR=${COMPARISON}_SNPs
mkdir ./$OUT_DIR


### GET FST
# To get fsts of multiple populations against each other (ie, in a table of fsts), you need to run them pairwise.
# The output file ("${COMPARISON}.weir.fst") will have the fst of each SNP.
echo "Finding pairwise Weir-Cockerham Fst with vcftools."
vcftools --vcf $VCF_FILE --weir-fst-pop $POP_FILE1 --weir-fst-pop $POP_FILE2 --out $OUT_DIR/$COMPARISON

## finding pairwise fixed SNPs for ID panel ##
#Vgrep to remove all rows with "-nan", negative fst, and decimals
echo "Preprocessing raw weir.fst output"
grep -v '-' $OUT_DIR/${COMPARISON}.weir.fst > $OUT_DIR/${COMPARISON}_SNPs.txt #Removes lines without an Fst value
grep -v $'\t0' $OUT_DIR/${COMPARISON}_SNPs.txt > $OUT_DIR/${COMPARISON}_fixed-SNPs.txt #Removes all Fst values that are not completely fixed

echo "Outputting all fixed SNPs"
cp $OUT_DIR/${COMPARISON}_fixed-SNPs.txt $OUT_DIR/${COMPARISON}_high-effect-SNPs.txt

echo "Outputting all SNPs with fixation > $MIN_FST_VALUE"
grep $'\t0.9' $OUT_DIR/${COMPARISON}_SNPs.txt >> $OUT_DIR/${COMPARISON}_high-effect-SNPs.txt
grep $'\t0.8' $OUT_DIR/${COMPARISON}_SNPs.txt >> $OUT_DIR/${COMPARISON}_high-effect-SNPs.txt
grep $'\t0.7' $OUT_DIR/${COMPARISON}_SNPs.txt >> $OUT_DIR/${COMPARISON}_high-effect-SNPs.txt

echo "Script is terminating."
