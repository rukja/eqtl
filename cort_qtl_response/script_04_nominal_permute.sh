#!/bin/bash
#$ -cwd


set -e

pwd
. /u/local/Modules/default/init/modules.sh


module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load R
module load plink

# Assign the task id to chromosome_number
if [ "$SGE_TASK_ID" -eq 23 ]; then
    chromosome_number="X"
else
    chromosome_number="$SGE_TASK_ID"
fi


echo "Chromosome number is: $chromosome_number"

condition=$1

vcf_file=$2

if [ -n "$condition" ]; then
    echo "Running analysis with condition: $condition"
else
    echo "Running control analysis"
fi

directory="./chr$chromosome_number"

mkdir -p "$directory"

chr=${chromosome_number}
if [ -z "$condition" ]; then
                my_bed="$directory/ctrl_chr${chromosome_number}.bed"
                filtered_name="filtered_nominals_${chr}.txt"
                variant_name="variant_list_${chr}.txt"
                plink_outputs="plink_outputs"
                nominals_name="$directory/nominals_${chromosome_number}_1.txt"
                permute_name="$directory/variants_permute_${chromosome_number}.txt"
                output_prefix="$directory/fdr_${chromosome_number}"
                conditional_name="$directory/conditional_pass_${chromosome_number}.txt"
                permutation_dir="ctrl_permutation_dir"
                cov_file="pca_bed_output_ctrl/ctrl_optimal_pc_cov.txt"
									
		zcat ./bigBed/ctrl.bed.gz | head -n 1 > "$my_bed" # Extract header
		tabix ./bigBed/ctrl.bed.gz "$chromosome_number" >> "$my_bed"
else
                my_bed="$directory/${condition}_chr${chromosome_number}.bed"
            	filtered_name="${condition}_filtered_nominals_${chr}.txt"
                variant_name="${condition}_variant_list_${chr}.txt"
                plink_outputs="${condition}_plink_outputs"
                nominals_name="$directory/${condition}_nominals_${chromosome_number}_1.txt"
                permute_name="$directory/${condition}_variants_permute_${chromosome_number}.txt"
                output_prefix="$directory/${condition}_fdr_${chromosome_number}"
                conditional_name="$directory/${condition}_conditional_pass_${chromosome_number}.txt"
                permutation_dir="${condition}_permutation_dir"
                cov_file="pca_bed_output_${condition}/${condition}_optimal_pc_cov.txt"

		zcat ./bigBed/${condition}.bed.gz | head -n 1 > "$my_bed" # Extract header
		tabix ./bigBed/${condition}.bed.gz "$chromosome_number" >> "$my_bed"
fi

if [ "$chromosome_number" = "X" ]; then
    sed -i 's/^X\t/23\t/' "$my_bed"
fi
# Extract chromosome level phenotype data

bgzip "$my_bed"
tabix -p bed "${my_bed}.gz"

# Create permutation folder

if [ ! -d ${permutation_dir} ]; then
  mkdir ${permutation_dir}
fi


# Extract chromosome level genotype data

my_vcf="$directory/chr${chromosome_number}.vcf"


# Check if sample_names.txt exists and is not empty
if [ ! -s sample_names.txt ]; then
    echo "Error: sample_names.txt is missing or empty"
    exit 1
fi

# Extract data with error checking
if ! tabix -h "${vcf_file}" "$SGE_TASK_ID" > "$my_vcf"; then
    echo "Error: tabix extraction failed"
    exit 1
fi

# Verify the VCF was created
if [ ! -s "$my_vcf" ]; then
    echo "Error: Output VCF is empty or was not created"
    exit 1
fi

if ! bgzip "$my_vcf"; then
    echo "Error: bgzip compression failed"
    exit 1
fi

if ! tabix -p vcf "${my_vcf}.gz"; then
    echo "Error: tabix indexing failed"
    exit 1
fi




# Perform nominal pass analysis (EDIT THE THRESHOLD AFTER --nominal)

QTLtools cis --vcf "${my_vcf}.gz" --bed "${my_bed}.gz" --cov ${cov_file} --nominal 1 --window 10000 --std-err --out "$directory/header_${chromosome_number}.txt" --chunk 0 1

QTLtools cis --vcf "${my_vcf}.gz" --bed "${my_bed}.gz" --cov ${cov_file} --nominal 1 --window 10000 --std-err --out "$directory/nominals_${chromosome_number}.txt" --chunk 1 1

cat "$directory/header_${chromosome_number}.txt" "$directory/nominals_${chromosome_number}.txt"  > "${nominals_name}"
rm "$directory/header_${chromosome_number}.txt"
rm "$directory/nominals_${chromosome_number}.txt"


# Perform permutation pass analysis

QTLtools cis --vcf "${my_vcf}.gz" --bed "${my_bed}.gz" --cov ${cov_file} --permute 5000 --window 10000 --std-err --out "$directory/permute_header_${chromosome_number}.txt" --chunk 0 1

QTLtools cis --vcf "${my_vcf}.gz" --bed "${my_bed}.gz" --cov ${cov_file} --permute 5000 --window 10000 --std-err --out "$directory/permute_${chromosome_number}.txt" --chunk 1 1 --seed 19

cat "$directory/permute_header_${chromosome_number}.txt" "$directory/permute_${chromosome_number}.txt" > "${permute_name}"


rm "$directory/permute_header_${chromosome_number}.txt"
rm "$directory/permute_${chromosome_number}.txt"

cp ${permute_name} ${permutation_dir}

echo "----------------------------"
echo "Done"

