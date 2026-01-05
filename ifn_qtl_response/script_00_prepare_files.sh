#!/bin/bash
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu
#$ -cwd
#$ -l h_rt=1:00:00,h_data=2G
#$ -pe shared 1
#$ -o joblog.$JOB_ID
#$ -j y

set -e

pwd
. /u/local/Modules/default/init/modules.sh


module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12


cd bigBed
bgzip *


touch ../beds.txt
ls >> ../beds.txt



tabix -p bed *

cd ..
zcat $(ls ./bigBed/*.bed.gz | head -n1) | head -n1 | cut -f7- | tr '\t' '\n' | sort | uniq > ./bigVCF/sample_names.txt
zcat $(ls ./bigBed/*.bed.gz | head -n1) | head -n1 | cut -f7- | tr '\t' '\n' | sort | uniq > ./sample_names.txt


echo "BED preparation complete."
