!/bin/bash
###################################################################################################
# This is a WIP script designed to run PCAs from the output files of 'find_pairwise_snps_high_fst.sh'
# Much of the source code is taken from 'KMJ-fixed-SNPs.txt'
# OUTPUT: <project_name.eigenvec> <project_name.eigenval>
###################################################################################################
### FIXED SNPS POST "FIND-PAIRWISE_SNPS_HIGH_FST.sh"
### This block of text greps out the header from each set of high-effect-SNPs

DATASET=Houbara_pairwise_sort
for SNPfile in ./*/*high-effect-SNPs.txt; do
PREFIX=`basename $SNPfile | sed 's/.txt//g` #This does not work as designed. 
grep -v 'CHROM' $SNPfile | sort >> ${DATASET}_sort.txt
done
HIEFF_SNPFILE=${DATASET}_sort.txt

# This line pulls the high-effect SNPs out of the VCF
vcftools --vcf $VCF_FILE --positions $HIEFF_SNPFILE --recode --recode-INFO-all --out ${DATASET}_high-effect-SNPs
HIEFF_VCF=${DATASET}_high-effect-SNPs.recode.vcf
HIEFF_PREFIX=`echo $HIEFF_VCF | sed 's/.recode.vcf//g'`

# This line converts the high-effect SNPs to a bedfile
plink --vcf $HIEFF_VCF --allow-extra-chr --double-id --make-bed --out $HIEFF_PREFIX
# This generates an eigenvalue and eigenvector file, for analysis in rstudio.
plink --bfile $HIEFF_PREFIX --pca --allow-extra-chr --out $HIEFF_PREFIX

setwd("/Users/kcollier/Downloads/")
library("ggplot2")
pca<- read.table("{$HIEFF_PREFIX}.eigenvec", sep=" ", header=F)

ggplot(data=pca, aes(x=V3, y=V4, color=V2)) + geom_point() + geom_text(size=2,aes(label=V1, size=NULL), nudge_y=.02) + theme(legend.position="none")
#These are saved in "PCAs/EUROTRIP3"
