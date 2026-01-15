#!/bin/bash
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu
#$ -cwd
#$ -l h_rt=2:00:00,h_data=2G
#$ -pe shared 2
#$ -o joblog_pc_bed.$JOB_ID
#$ -j y

set -e

# Check for input file parameter
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 input_file.txt number_of_genotype_PCs path_to_vcf covariate_inclusion"
    echo "Where input_file.txt contains one .bed.gz filename per line, and covariate_inclusion refers to whether the custom covariate fill should be included (Y or N)"
    exit 1
fi

INPUT_FILE="$1"

# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

vcf_pcs="$2"


# Check if vcf pcs are valid integers

re='^[0-9]+$'
if ! [[ $vcf_pcs =~ $re ]] ; then
   echo "Error: Provided VCF PCs, ${vcf_pcs} is not an integer"
   exit 1
fi

my_vcf="$3"
inclusion="$4"


# Check if vcf file exists
if [ ! -f "$my_vcf" ]; then
    echo "Error: VCF file '$my_vcf' not found"
    exit 1
fi

# Check if inclusion is valid
if [ "${inclusion}" != "Y" ] && [ "${inclusion}" != "N" ]; then
    echo "Error: information about covariate inclusion not present"
    exit 1
fi

pwd
. /u/local/Modules/default/init/modules.sh

module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load R
module load plink
module load vcftools/0.1.16

vcf_directory="./bigVCF"
bed_directory="./bigBed"




# Preparing for R script execution
module load apptainer
if [ ! -d "./hoffman_log_Rscripts" ]; then
  mkdir ./hoffman_log_Rscripts
fi
export R_LIBS_USER=$HOME/R/APPTAINER/h2-rstudio_4.4.0

bcftools query -l ${my_vcf} > vcf_samples.txt


# Process each file in the input list
while IFS= read -r bed_file; do
    echo "Processing $bed_file"

    # Extract original_bed name without .gz extension
    original_bed=$(basename "$bed_file" .bed.gz)

    echo "Original bed file: $original_bed"

    # Create output directory for this bed file
    output_directory="./pca_bed_output_${original_bed}"

    if [ ! -d "${output_directory}" ]; then
        mkdir "${output_directory}"
        echo "Created directory: ${output_directory}"
    fi

    my_bed="${bed_directory}/${bed_file}"

    # Obtain PCA results
    echo "Calculating PCs..."
    /usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.4.0.sif Rscript ./calculate_pcs.R \
        "${my_bed}" \
        "${output_directory}/genes.corrected" \
        > ./hoffman_log_Rscripts/find_pc_bed_log.$JOB_ID 2>&1

    # Create covariate file
    (head -n 1 "${vcf_directory}/PCs.pca"  && tail -n +2 "${vcf_directory}/PCs.pca"  | head -n ${vcf_pcs}) > "${vcf_directory}/PCs_selected.pca"
    cov_file="${output_directory}/${original_bed}.covariates.pcs.txt"

    echo "Covariate merging..."
    /usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.4.0.sif Rscript ./merge_pcs.R \
        "${vcf_directory}/PCs_selected.pca"  \
        "./covariates.txt" \
        "${output_directory}/genes.corrected.pca" \
        "${inclusion}" \
        ${cov_file} \
        > ./hoffman_log_Rscripts/merge_pcs_output.$JOB_ID 2>&1

    echo "Covariate investigation..."
    /usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.4.0.sif Rscript ./redundant_pcs.R \
        "$cov_file" \
        "${cov_file}" \
        > ./hoffman_log_Rscripts/cov_investigation_output.$JOB_ID 2>&1

    cp ${my_bed} "${bed_directory}/copy_${bed_file}"

    echo "Prepare order..."
    /usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.4.0.sif Rscript ./reorder_columns.R \
        "${my_bed}" \
        "${bed_directory}/${original_bed}.bed" \
        > ./hoffman_log_Rscripts/reorder_bed_pcs.$JOB_ID 2>&1

    sed -i '/^#/!s/^chr//' ${bed_directory}/${original_bed}.bed
	
    echo "Re-index BED"
    bgzip -f "${bed_directory}/${original_bed}.bed"
    tabix -p bed -f ${my_bed}

done < "$INPUT_FILE"



