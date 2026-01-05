#!/bin/bash
#$ -cwd
#$ -o extract_genotypes_variant_list.out.log
#$ -j y
#$ -l h_rt=6:00:00,h_data=2G
#  PLEASE CHANGE THE NUMBER OF CORES REQUESTED AS NEEDED:
#$ -pe shared 3
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu

set -e

pwd
. /u/local/Modules/default/init/modules.sh

module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load vcftools

# List of variants

variants_file="$1"

# Validate input files exist
if [ ! -f "$variants_file" ]; then
    echo "Variants file not found: $variants_file"
    exit 1
fi

header_file="header.txt"  

# Create an empty file to combine all genotypes 
output_file="./combined_genotypes_of_interest.tsv"
> $output_file  # Clear the output file if it exists

# Extract the header from the VCF and save it
zgrep "^#CHROM" "./chr1/chr1.vcf.gz" | cut -f1-5,10- | sed 's/^#//' > $header_file

echo "Header created"



while IFS=' ' read -r variant chr; do

	if [ -z "$chr" ]; then
		echo "Error: the input format must be variant chr. Chr not supplied!"
		exit 1
	fi

	if [ -z "$variant" ]; then
		echo "Error: the input format must be variant chr. Variant not supplied!"
		exit 1
	fi

	dir="./$chr"

	if [ ! -d "$dir" ]; then
		echo "Error: The directory does not exist. Please check the chromosome input. The form of input should be chr1, chr9, chrX."
		exit 1
	fi
	
	echo "Processing variant $variant in chromosome $chr"
	
	cd "$dir"

	echo "Extracting variant"
	vcftools --gzvcf "./${chr}.vcf.gz" --snp "$variant" --recode --out ${variant}

   	 # Use bcftools to query the VCF and get the genotypes
   	 bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT[\t%GT]\n' ${variant}.recode.vcf > ${variant}_genotypes.tsv

		if [ ! -s "${variant}_genotypes.tsv" ]; then
			echo "Error: error extracting genotypes"
			exit 1
		fi

		cat "${variant}_genotypes.tsv" >> "../${output_file}"

   	 # Remove the temporary VCF and TSV files after processing
   	 rm ${variant}.recode.vcf ${variant}_genotypes.tsv
	 cd ..
done < "$variants_file"

# Combine the header with the genotypes
cat $header_file $output_file > "combined_filtered_variants_genotypes.tsv"

# Clean up
rm $header_file $output_file

