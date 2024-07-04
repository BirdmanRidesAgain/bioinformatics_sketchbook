#!/bin/bash
# $1 is the prefix of the output file

# get all vcfs as elements in an array
vcfs=($(ls *vcf.gz))
DEPTH_OUTPUT="coverage_results.txt"

# create output file to show filters
echo -e "file\tprefix\tmean\t2.5th_Percentile\t97.5th_Percentile" > ${1}_${DEPTH_OUTPUT}

# Filtering loop:
for ((i=0; i<${#vcfs[@]}; i++)); do
    a=${vcfs[i]%.*}
    prefix=${a%.*} # gets us the prefix
    
    # creates array of first ~100000 lines of the vcf. could be improved by using whatever bcftools/samtools builtin probably exists
    ind_cov=($(zcat ${vcfs[i]} | head -100000 | grep -v '^#' | awk '{print $10}' FS='\t' OFS='\n' | cut -d: -f3 | grep -oE '[0-9]+'))
    ind_cov_len=${#ind_cov[@]} # gets the number of variants where we looked at depth.
    n=$ind_cov_len #copy ind_cov_len to value we can decrement
    sum=0 # initializes our accumulator

    # inner loop to get sum of all depths to calc mean:
    while [ $n -gt 0 ] # while n is greater than 0, run the loop.
    do 
        num=${ind_cov[`expr $n - 1`]} # get value from ind_cov
        sum=`expr $sum + $num` # add all values of ind_cov together
        n=`expr $n - 1` # decrement our counter
    done
    mean=`expr $sum / $ind_cov_len | bc -l`
    #fixme - these probably need to be calculated based on mean and not whatever Thierry is doing
    percent_2_5=$(echo "${ind_cov[*]}" | tr ' ' '\n' | sort -n | awk 'BEGIN{q=0.025} {a[NR]=$1} END {print a[int(NR*q)];}') # it rounds down
    percent_97_5=$( echo "${ind_cov[*]}" | tr ' ' '\n' | sort -n | awk 'BEGIN{q=0.975} {a[NR]=$1} END {print a[int(NR*q)];}')
    
    echo -e "${vcfs[i]}\t$prefix\t$mean\t$percent_2_5\t$percent_97_5" # so the user can see
    echo -e "${vcfs[i]}\t$prefix\t$mean\t$percent_2_5\t$percent_97_5" >> ${DEPTH_OUTPUT}

    # Add an MQ filter - for falcons 55
    # remove low-quality contigs    
    # Bcftools view - filters out multiallelic, missing, indel, low-GQ and abnormal depth data
    bcftools view -M2 -m2 -U -V indels ${vcfs[i]} | bcftools filter -i "FORMAT/GQ>30 && FORMAT/DP<${percent_97_5} && FORMAT/DP>${percent_2_5}" -o ${prefix}_filt.vcf.gz -Oz
done


# Merge all vcfs and filter on MAF and QUAL
for i in *_filt.vcf.gz; do
bcftools index $i 
done
ls *_filt.vcf.gz | xargs bcftools merge | bcftools filter -i "MAF>0.03" | bcftools view -o ${1}_filt.vcf.gz -Oz 
bcftools index ${1}_filt.vcf.gz

# print useful information for a human:
num_snps=`bcftools query -f "\n" ${1}_filt.vcf.gz | wc -l`
echo -e "There are $num_snps snps in ${1}_filt.vcf.gz"
