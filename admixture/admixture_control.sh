#!/bin/bash

USAGE="Usage: admixture_control.sh <input.bed> <input.fam> <num_iterations> <max_k>"
OUTPUTS="Outputs: all admixture output files to a directory called input_admixture/K<1-max_k>"
NOTES="
\tk-max should be your number of populations.
\tYou are assumed to have edited your .fam/.pop file beforehand. It should consist of:
Pop_1 Ind_1
Pop_2 Ind_2
Pop_2 Ind_3
Pop_3 Ind_4

All pop/ind pairs need to be on separate lines, separated by whitespace.
"
if [ -z $1 ]
then
    echo -e "$USAGE\n${OUTPUTS}"
    echo -e "$NOTES"
    exit 1
fi

####################################################################################
### SOURCING VARIABLES AND CONDA
eval "$(conda shell.bash hook)"
BED=$1
PREFIX=$(basename -s .bed "${BED}")
POP=$2

ITERATIONS=$3
MAX_K=$4
BS_REPS=10

##################################
# CREATE TOP-LEVEL DIRECTORY
ADMIXTURE_OUTDIR=./${PREFIX}_admixture
if [ -d "$ADMIXTURE_OUTDIR" ]; then
    rm -r $ADMIXTURE_OUTDIR
fi
mkdir ${ADMIXTURE_OUTDIR}

##################################
# RUN ADMIXTURE
conda activate admixture

# We are overwriting the ${PREFIX}.fam file admixture attempts to find b/c the overwritten file will never have pop. information.
# The biologist needs to explicitly provide that info.
cp $POP ${PREFIX}.fam
echo "Beginning admixture runs for K=1 to K=$MAX_K"
admixture -B${BS_REPS} -j${THREADS} --seed=$j $BED $i #> ${PREFIX}.K${i}.rep${j}.log

for i in $(seq 1 $MAX_K); do
    echo "  Beginning K=$i"
    mkdir -p $ADMIXTURE_OUTDIR/K${i}
    K_DIR=$ADMIXTURE_OUTDIR/K${i}

    # repeat 30x, each with a different seed
    for j in $(seq 1 $NUM_ITERATIONS); do
        echo "      K=$i, replicate $j"
        admixture -B${BS_REPS} -j${THREADS} --seed=$j $BED $i > ${PREFIX}.K${i}.rep${j}.log

        # move output files to the appropriate directory
        LOG=${PREFIX}.K${i}.rep${j}.log
        PQ_PREFIX=${PREFIX}.K${i}.rep${j}

        mv $LOG $K_DIR
        mv ${PREFIX}.${i}*.P $K_DIR/${PQ_PREFIX}.P
        mv ${PREFIX}.${i}*.Q $K_DIR/${PQ_PREFIX}.Q
        mv ${PREFIX}.${i}*.bias $K_DIR/${PQ_PREFIX}.bias
        mv ${PREFIX}.${i}*.se $K_DIR/${PQ_PREFIX}.se
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
        CURRENT_LOG=$ADMIXTURE_OUTDIR/K${i}/${PREFIX}.K${i}.rep${j}.log
        LOGLIKE=`cat $CURRENT_LOG | grep -v "(delta)" | grep Loglikelihood | awk '{print $NF}'`
        echo $i $LOGLIKE >> $ADMIXTURE_OUTDIR/loglike.K${i}.txt
    done
    echo "  Finished collecting log likelihood values for K=$i"
done

SUMMARY_FILE=$ADMIXTURE_OUTDIR/loglike.K1-K${MAX_K}.txt
echo "num_K" "log_likelihood" > $SUMMARY_FILE
cat $ADMIXTURE_OUTDIR/loglike.K*.txt >> $SUMMARY_FILE
echo "Finished collecting log likelihood values for K=1 to K=$MAX_K"
