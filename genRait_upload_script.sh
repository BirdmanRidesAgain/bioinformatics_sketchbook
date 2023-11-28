#!/bin/bash
# Usage: 'bash genRait_upload_script.sh <file_or_directory_to_be_uploaded> <project_name>'
WDIR=`pwd`;
GEN=$WDIR/gen-cli;
UPLOAD=$1;
PROJECT_NAME=$2;

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
echo "$p" > job.txt;


$JOB=${PROJECT_NAME}_upload_$(date +%s).sh;
cat job.txt >> $JOB;
chmod +x $JOB;
rm job.txt;
sh $JOB;
