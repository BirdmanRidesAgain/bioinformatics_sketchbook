#/bin/bash
# Houbara bash commands
# ---------------------
# This script is used to run bcftools variantcaller on a subset of Houbara
# The refseq is houbaraONT_assembly-racon-haplotigs.fa
WD=/mnt/data-2/Houbara
REFSEQ=${WD}/refseq/houbaraONT_assembly-racon-haplotigs.fa
ARRAY_NOVCF=("H4348" "M03N03872" "M03N03943" "M05N06139" "M12N18705" "N15390" "N2766" "N2836" "N3469" "N4088" "N4090" "N4098" "N4391" "P2595" "P2630" "P2945" "P2946" "P3796" "U87M00087")

for i in "${ARRAY_NOVCF[@]}"; do
    echo $i + " is being processed"
    bcftools mpileup -Ou -f $REFSEQ \
    -a FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \
     ${WD}/raw-reads/${i}/${i}_CL-SO-RG-MD.bam | \
     bcftools call -mv -Oz \
     -f GQ,GP \
     -o ${WD}/vcfs/${i}.vcf.gz
done

