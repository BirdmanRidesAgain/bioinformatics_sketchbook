#!/bin/bash
# usage: psmc_module.sh <assembly.fq.gz> <bs_reps> <project_name>
# outputs:

########################################################################
### SETTING VARIABLES
INPUT_FQ=$1
BS_REPS=$2
PROJECT_NAME=$3

########################################################################
### PREPARING INPUT FILE
# Get the input file with 'input_file_prep.sh'

########################################################################
### PSMC CALCULATION
### Without bootstrapping


### With bootstrapping
fq2psmcfa -q 20 $INPUT_FQ > ${PROJECT_NAME}.psmcfa     ## format assembly to psmcfa format
splitfa ${PROJECT_NAME}.psmcfa > ${PROJECT_NAME}_split.psmcfa       ## split long chromosome sequences to short ones
psmc -o ${PROJECT_NAME}.psmc ${PROJECT_NAME}.psmcfa     ## output psmc results from un-split psmcfa file
seq ${BS_REPS} | xargs -i -P 10 echo psmc -o ${PROJECT_NAME}_round-{}.psmc ${PROJECT_NAME}_split.psmcfa | sh        ## use "-P <int>" to set maximum # of processes run at once. -p for interactive mode.

cat ${PROJECT_NAME}.psmc round-*.psmc ${PROJECT_NAME}_split.psmcfa > ${PROJECT_NAME}_${BS_REPS}-BSreps.psmc
psmc_plot.pl -p -Y 50000 ${PROJECT_NAME}_${BS_REPS}-BSreps_out
