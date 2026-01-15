#!/bin/bash
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu
#$ -cwd
#$ -l h_rt=23:00:00,h_data=8G
#$ -pe shared 4
#$ -o joblog.$JOB_ID
#$ -j y

set -e

# Check for input file parameter
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input_file.txt vcf_name.hwe.vcf.gz vcf_pcs"
    echo "Where input_file.txt contains one .bed.gz filename per line, vcf_name.hwe.vcf.gz is the path to the VCF, and vcf_pcs contains the number of PCs to correct in the VCF"
    exit 1
fi

INPUT_FILE="$1"
my_vcf="$2"


# Check if input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found"
    exit 1
fi

# Check if VCF file exists
if [ ! -f "$my_vcf" ]; then
    echo "Error: VCF file '$my_vcf' not found"
    exit 1
fi

# Check if vcf pcs are valid integers
vcf_pcs="$3"
re='^[0-9]+$'
if ! [[ $vcf_pcs =~ $re ]] ; then
   echo "Error: Provided VCF PCs, ${vcf_pcs} is not an integer"
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

    echo "PCs,Significant_eQTLs" > ${output_directory}/pc_optimization_results.csv

    cov_file="${output_directory}/${original_bed}.covariates.pcs.txt"

    total_lines=$(wc -l < "$cov_file")
    num_pcs_available=$((total_lines - 1))
    # Iterate and perform analyses
    for num_pcs in $(seq ${vcf_pcs} ${num_pcs_available}); do
        num_other_pcs=$((num_pcs - vcf_pcs))
        echo "Testing with ${vcf_pcs} VCF PCs + ${num_other_pcs} bed PCs (${num_pcs} total)..."

        # Create a covariates file with only the first N additional PCs
        head -n 1 ${cov_file} > ${output_directory}/genes.corrected.${num_pcs}pcs.pca
        awk -v n=${num_pcs} 'NR > 1 && NR <= n+1' ${cov_file} >> ${output_directory}/genes.corrected.${num_pcs}pcs.pca

        # Run QTL mapping
        QTLtools cis --vcf ${my_vcf} \
                    --bed ${my_bed} \
                    --out ${output_directory}/qtl_results_${num_pcs}pcs.txt \
                    --nominal 0.05 \
                    --cov ${output_directory}/genes.corrected.${num_pcs}pcs.pca \
                    --window 10000


        # Count number of significant eQTLs
        significant_hits=$(awk '$12 < 0.05' ${output_directory}/qtl_results_${num_pcs}pcs.txt | wc -l)
        echo "Number of significant eQTLs with ${num_pcs} PCs: ${significant_hits}"
        echo "${num_pcs},${significant_hits}" >> ${output_directory}/pc_optimization_results.csv
    done

    # Find optimal number of PCs
    optimal_pcs=$(sort -t, -k2,2nr ${output_directory}/pc_optimization_results.csv | head -n1 | cut -d, -f1)
    echo "Optimal number of PCs for ${original_bed}: ${optimal_pcs}"

    bed_pc="${output_directory}/genes.corrected.${optimal_pcs}pcs.pca"

    mv "${output_directory}/qtl_results_${optimal_pcs}pcs.txt" "${output_directory}/${original_bed}_optimal_results.txt"

    mv $bed_pc "${output_directory}/${original_bed}_optimal_pc_cov.txt"

    echo "Completed processing ${bed_file}"
    echo "----------------------------------------"

done < "$INPUT_FILE"

cd ${bed_directory}
#rm temp*

cd ..

echo "All files processed successfully"
