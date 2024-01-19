#!/bin/bash
# Usage:  ./admixture_control_script.sh <input.bed> <num_iterations> <max_k>
# Outputs: all admixture output files to a directory called <prefix>_output/K<1-max_k>

##################################
# ACTIVATE CONDA
conda activate admixture

##################################
### SETTING VARIABLES
#User-supplied vars
INPUT_BED=$1
NUM_ITERATIONS=$2
MAX_K=$3

#Other vars
WD=`pwd`
THREADS=40
PREFIX=`ls $INPUT_BED | sed 's/.bed//g'`

##################################
# CREATE DIRECTORIES
ADMIXTURE_OUTDIR=$WD/${PREFIX}_output

##################################
# RUN ADMIXTURE
# Loop through each value of K
echo "Beginning admixture runs for K=1 to K=$MAX_K"

for i in $(seq 1 $MAX_K); do
    echo "  Beginning K=$i"
    mkdir -p $ADMIXTURE_OUTDIR/K${i}
    K_DIR=$ADMIXTURE_OUTDIR/K${i}

    # repeat 30x, each with a different seed
    for j in $(seq 1 $NUM_ITERATIONS); do
        echo "      K=$i, replicate $j"
        admixture --cv $INPUT_BED -j${THREADS} --seed=$j $i > ${PREFIX}.K${i}.rep${j}.log

        # move output files to the appropriate directory
        LOG=${PREFIX}.K${i}.rep${j}.log
        PQ_PREFIX=${PREFIX}.K${i}.rep${j}

        mv $LOG $K_DIR
        mv ${PREFIX}.${i}*.P $K_DIR/${PQ_PREFIX}.P
        mv ${PREFIX}.${i}*.Q $K_DIR/${PQ_PREFIX}.Q
    done
    echo "  Finished K=$i"
done
echo "Finished admixture runs for K=1 to K=$MAX_K"


##################################
# GET LOG LIKELIHOOD VALUES FROM ALL VALUES OF K
# Loop through each value of K
echo "Collecting log likelihood values for K=1 to K=$MAX_K"
for i in $(seq 1 $MAX_K); do
    echo "  Collecting log likelihood values for K=$i"
    # find the log likelihood value for each replicate
    for j in $(seq 1 $NUM_ITERATIONS); do
        echo "      K=$i, replicate $j"
        LOGLIKE=`cat $ADMIXTURE_OUTDIR/K${i}/${PREFIX}.K${i}.rep${j}.log | grep -v "(delta)" | grep Loglikelihood | awk '{print $NF}'`
        echo $LOGLIKE >> $ADMIXTURE_OUTDIR/loglike.K${i}.txt
    done
    echo "  Finished collecting log likelihood values for K=$i"
done
echo "Finished collecting log likelihood values for K=1 to K=$MAX_K"


##################################
# FINDING BEST K WITH EVANNO METHOD
    # do stuff here.


##################################
# TERMINATE PROGRAM
echo "Finished admixture_control_script.sh"



##################################
# SAMPLE CODE FOR LATER SPLITS - THIERRY, PLEASE IGNORE
### SPLIT ONE ###

# command to convert fam into an ind_list
#cat *fam | grep -v "Yemen" | sed 's/ 0 0 0 -9//g' > ind_list_first_split.txt

# command to use ind_list file to remove the population from the bed file
#plink --bfile HoubaraFeb23_noUndulata_hiqual_admixture --allow-extra-chr --chr-set 95 --keep ind_list_first_split.txt --make-bed --out HoubaraFeb23_noUndulata_first_split

### SPLIT TWO ###
# Repeat the above steps for the second split


#mv $i $(echo $i | sed 's/Undulate/Undulata/g')

# Create directories K1 - K9 in bash
#for i in {1..9}; do mkdir K$i; done



###################
#houbaraset=allhoubara
#num_iteration=30
#unique seed needed (1:30)
#range_k=1:9

# For each of the 270 runs, extract the converged loglikehihood value
    # For each value of k, write all 30 loglikihood values to a file
        # get average of the 30 loglike values

# Apply Evanno method to the output file to get the optimal k
    # Find the delta K value for each k
    # Graph two lines:

