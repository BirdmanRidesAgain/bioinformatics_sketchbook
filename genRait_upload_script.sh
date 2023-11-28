#!/bin/bash
# Usage: 'bash genRait_upload_script.sh <file_or_directory_to_be_uploaded> <project_name>'

########################################################################
### SETTING VARIABLES
WDIR=`pwd`;
UPLOAD=$1;
PROJECT_NAME=$2;
GEN=$WDIR/gen-cli;

########################################################################
### CONFIGURING GENRAIT
echo "Configuring GenRait login";
$GEN config kcollier@reneco.org B\!rdman2023;
echo "Generating metadata";
$GEN metadata $UPLOAD -r -c 100;

NL=$'\n';
ALLRUNS=$(ls $WDIR/runs);
 
for RUN in $ALLRUNS
do
p+="$GEN upload -m $WDIR/runs/$RUN -h -b 32 -p /${PROJECT_NAME}/ ;sleep 2;"$'\n'
done 
echo "$p" > ${PROJECT_NAME}.txt;


JOB=${PROJECT_NAME}_upload_$(date +%s).sh;
cat ${PROJECT_NAME}.txt >> $JOB;
chmod +x $JOB;
rm ${PROJECT_NAME}.txt;
sh $JOB;
