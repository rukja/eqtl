#!/bin/bash
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu
#$ -cwd
#$ -l h_rt=2:00:00,h_data=20G
#$ -o joblog.$JOB_ID
#$ -j y


set -e

pwd
. /u/local/Modules/default/init/modules.sh

module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load R
module load plink
module load vcftools/0.1.16

## Required inputs ##
### a path to the current vcf: my_vcf
### a minor allele frequency: MAF



# Check if required arguments are provided
if [ $# -lt 2 ]; then
    echo "Usage: $0 <vcf_path> <MAF>"
    echo "  <vcf_path>: Path to the input VCF file"
    echo "  <MAF>: Minor allele frequency threshold (e.g., 0.01)"
    exit 1
fi

# Assign command line arguments to variables
my_vcf="$1"
MAF="$2"

## Begin analysis

vcf_base=$(basename "$my_vcf" | sed -E 's/\.vcf(\.bgz)?$//')


temp_vcf="${vcf_base}.temp.vcf"
filtered_vcf="${vcf_base}.filtered.vcf"
hwe_vcf="${vcf_base}.hwe"
prune_output="${vcf_base}.pruned"


echo "Starting preparation of input files"
echo "Using VCF: $my_vcf"
echo "Using MAF threshold: $MAF"

# Check if VCF file exists
if [ ! -s "$my_vcf" ]; then
    echo "Error: Input VCF file '$my_vcf' is missing or empty"
    exit 1
fi

if [ ! -s sample_names.txt ]; then
    echo "Error: sample_names.txt is missing or empty"
    exit 1
fi

variant_count_before=$(bcftools view -S sample_names.txt "$my_vcf" | \
                       bcftools view -H | wc -l)
echo "Variants before genotype filtering: $variant_count_before"

# VCF subset based on patients with gene expression and MAF frequency
temp_filled="${vcf_base}.filled.vcf.gz"

# Step 1: Subset samples and apply MAF filter
if ! bcftools view -S sample_names.txt "$my_vcf" | \
     bcftools view -i "MAF[0]>=$MAF" | \
     bcftools +fill-tags -Oz -o "$temp_filled" -- -t AC_Hom,AC_Het,AN; then
    echo "Error: bcftools sample selection and fill-tags failed"
    exit 1
fi

# Index the temporary file
tabix -p vcf "$temp_filled"

# Step 2: Apply genotype count filters (all individuals >= 3)
if ! bcftools view -i 'AC_Hom>=6 && AC_Het>=3 && AN >= AC_Hom + 2*AC_Het + 6' \
     "$temp_filled" -o "$temp_vcf"; then
    echo "Error: bcftools genotype filtering failed"
    rm -f "$temp_filled" "${temp_filled}.tbi"
    exit 1
fi

# Clean up temporary file
rm -f "$temp_filled" "${temp_filled}.tbi"

variant_count_after=$(bcftools view -H "$temp_vcf" | wc -l)
echo "Variants after genotype filtering: $variant_count_after"
echo "Variants removed: $((variant_count_before - variant_count_after))"

mv "$temp_vcf" "$filtered_vcf"


bcftools query -f '%CHROM\t%POS\t%AC_Hom\t%AC_Het\t%AN\n' ${filtered_vcf} | head -20 > genotype_filter_check.txt


# VCF hwe correction
# Convert filtered VCF to PLINK binary files
plink --vcf "${filtered_vcf}" --const-fid 0 --make-bed --out "${filtered_vcf}_plink"

plink --bfile "${filtered_vcf}_plink" --missing --out pre_filter_check
echo "Individual missingness rates:"
cat pre_filter_check.imiss
echo "Variant missingness rates summary:"
cat pre_filter_check.lmiss | head -20

# Apply filtering in PLINK:
# - HWE p-value threshold 1e-6
# - Remove variants with >=5% missingness (--geno 0.05)
# - Filter out individuals with >5% missing SNPs
plink --bfile "${filtered_vcf}_plink" \
  --hwe 1e-6 \
  --geno 0.05 \
  --mind 0.10 \
  --mac 14 \
  --make-bed --out "${hwe_vcf}_plink_filtered"

  # Convert filtered PLINK files back to VCF
plink --bfile "${hwe_vcf}_plink_filtered" --recode vcf-iid --out "${hwe_vcf}"

# Compress the final VCF
bgzip -f "${hwe_vcf}.vcf"

# Check if the compressed VCF is not empty
if [ ! -s "${hwe_vcf}.vcf.gz" ]; then
    echo "Error: Output VCF file '${hwe_vcf}.vcf.gz' is empty after filtering"
    exit 1
fi

# Validate VCF format
if ! bcftools view -h "${hwe_vcf}.vcf.gz" &>/dev/null; then
    echo "Error: The filtered VCF file '${hwe_vcf}.vcf.gz' appears malformed"
    exit 1
fi

# Index the VCF
tabix -p vcf "${hwe_vcf}.vcf.gz"


# VCF LD pruning with PLINK for subsequent PCA calculation

echo "Starting PLINK pruning"

if [ ! -s exclusion_regions_hg19.txt ]; then
    echo "Error: exclusions file is missing or empty"
    exit 1
fi

plink --vcf "${hwe_vcf}.vcf.gz" --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions_hg19.txt
plink --vcf "${hwe_vcf}.vcf.gz" --extract plink.prune.in --make-bed --out ${prune_output}

# PCA calculation

echo "Starting PCA caluclation"

if [ ! -s flashpca.R ]; then
    echo "Error: R file is missing or empty"
    exit 1
fi

module load apptainer
if [ ! -d "./hoffman_log" ]; then
  mkdir ./hoffman_log
fi
export R_LIBS_USER=$HOME/R/APPTAINER/h2-rstudio_4.4.0
/usr/bin/time -v apptainer exec $H2_CONTAINER_LOC/h2-rstudio_4.4.0.sif Rscript ./flashpca.R "${prune_output}" > ./hoffman_log/output.$JOB_ID 2>&1

rm plink*
rm *pruned*
rm *filtered*
rm pre_*

echo "All analysis complete!"


