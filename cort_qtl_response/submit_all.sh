#!/bin/bash
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu
#$ -cwd
#$ -l h_rt=8:00:00,h_data=20G
#$ -pe shared 3
#  PLEASE CHANGE THE NUMBER OF CORES REQUESTED AS NEEDED:
#$ -t 1-23:1

#$ -o joblog.$JOB_ID.$TASK_ID
#$ -j y




# echo job info on joblog:
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID started on:   " `date `
echo " "

set -e

pwd
. /u/local/Modules/default/init/modules.sh

export R_LIBS_USER=$HOME/R/APPTAINER/h2-rstudio_4.4.0

my_vcf=$1

# Check if input file exists
if [ ! -f "$my_vcf" ]; then
    echo "Error: Input file '$my_vcf' not found"
    exit 1
fi

while read condition; do
    echo "Processing condition: $condition"
    /usr/bin/time -v ./script_04_nominal_permute.sh $condition ${my_vcf}
done < analyses.txt


echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `hostname -s`
echo "Job $JOB_ID.$SGE_TASK_ID ended on:   " `date `
echo " "


