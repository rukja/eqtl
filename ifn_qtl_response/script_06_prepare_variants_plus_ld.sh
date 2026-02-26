#!/bin/bash
#$ -cwd
#$ -o ld.out.log
#$ -j y
#$ -l h_rt=6:00:00,h_data=2G
#  PLEASE CHANGE THE NUMBER OF CORES REQUESTED AS NEEDED:
#$ -pe shared 3
# Notify when
#$ -m bea
#$ -M rithwiknarendra@ucla.edu

pwd
. /u/local/Modules/default/init/modules.sh

module load apptainer
module load qtltools/1.3.1
module load bcftools/1.11
module load htslib/1.12
module load R
module load plink

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <plink_prefix> <output_prefix>"
    exit 1
fi

PLINK_PRE=$1

OUTPUT_PRE=$2


## Process conditional hits and create files for extracting genotypes and independent hits

awk '!a[$0]++ {print}' ./chr*/chr*_conditionals.txt > conditionals_full.txt

echo "There are $(wc -l < conditionals_full.txt) total conditional hits"

awk 'NR==1 || $21==1' conditionals_full.txt > conditional_top_variants.txt

echo "There are $(wc -l < conditional_top_variants.txt) independent conditional hits"

awk 'NR>1 {print $8, "chr"$2}' conditionals_full.txt > variants.txt

awk 'NR>1 {print $1}' conditional_top_variants.txt | sort -u > reGenes.txt

echo "There are $(wc -l < reGenes.txt) reGenes"

## Perform LD analysis

echo "Running LD per gene..."

echo "gene,SNP_A,SNP_B,R2" > ${OUTPUT_PRE}_ALL_LD.csv

while read g; do

    echo "Processing gene: $g"

    # Extract SNP IDs (column 8) for this gene
    awk -v gene="$g" '$1==gene {print $8}' conditionals_full.txt > ${g}_snps.txt

    # Skip if fewer than 1 SNPs
    snp_count=$(wc -l < ${g}_snps.txt)
    if [ "$snp_count" -lt 2 ]; then
        echo "Skipping $g (only $snp_count SNP)"
        rm ${g}_snps.txt
        continue
    fi

    # Run LD
    plink \
        --bfile ./bigVCF/downstream_ld_output/${PLINK_PRE} \
        --extract ${g}_snps.txt \
        --r2 \
        --ld-window 99999 \
        --ld-window-kb 1000 \
        --ld-window-r2 0 \
        --out ${OUTPUT_PRE}_${g}_LD

    awk -v gene="$g" 'NR>1 {print gene","$3","$6","$7}' \
    ${OUTPUT_PRE}_${g}_LD.ld >> ${OUTPUT_PRE}_ALL_LD.csv

    rm ${OUTPUT_PRE}_${g}_LD.*

    rm ${g}_snps.txt

done < reGenes.txt



