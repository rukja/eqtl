#!/bin/bash
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu
#$ -cwd
#$ -l h_rt=1:00:00,h_data=8G
#$ -pe shared 4
#$ -o joblog.$JOB_ID
#$ -j y

set -e

pwd
. /u/local/Modules/default/init/modules.sh

module load apptainer
module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load R
module load plink

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 condition condition with all permutation results"
    exit 1
fi

PERMUTATION_CONDITION="$1"

head -n 1 -q ${PERMUTATION_CONDITION}_permutation_dir/*.txt | head -n 1 > ${PERMUTATION_CONDITION}_permutation_dir/permutations.txt
tail -n +2 -q ${PERMUTATION_CONDITION}_permutation_dir/*.txt >> ${PERMUTATION_CONDITION}_permutation_dir/permutations.txt


export R_LIBS_USER=$HOME/R/APPTAINER/h2-rstudio_4.4.0

/usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.4.0.sif Rscript ./qtltools_runFDR_cis.R \
    "${PERMUTATION_CONDITION}_permutation_dir/permutations.txt" \
    0.05 \
    "${PERMUTATION_CONDITION}_permutation_dir/thresholds" \
    > ./hoffman_log_Rscripts/output_thresholds.$JOB_ID



cov_file="pca_bed_output_${PERMUTATION_CONDITION}/${PERMUTATION_CONDITION}_optimal_pc_cov.txt"

for num in {1..22} X; do


    QTLtools cis --vcf "chr${num}/chr${num}.vcf.gz" --bed "chr${num}/${PERMUTATION_CONDITION}_chr${num}.bed.gz" --cov ${cov_file} --window 10000 --mapping "${PERMUTATION_CONDITION}_permutation_dir/thresholds.thresholds.txt" --out "chr${num}/header_conditional_pass_${num}.txt" --chunk 0 1
    QTLtools cis --vcf "chr${num}/chr${num}.vcf.gz" --bed "chr${num}/${PERMUTATION_CONDITION}_chr${num}.bed.gz" --cov ${cov_file} --window 10000 --mapping "${PERMUTATION_CONDITION}_permutation_dir/thresholds.thresholds.txt" --out "chr${num}/data_conditional_pass_${num}.txt" --chunk 1 1


    cat "chr${num}/header_conditional_pass_${num}.txt" "chr${num}/data_conditional_pass_${num}.txt" > "chr${num}/chr${num}_conditionals.txt"

    rm "chr${num}/header_conditional_pass_${num}.txt"
    rm "chr${num}/data_conditional_pass_${num}.txt"


    echo "Done with passes ${num}"


done






