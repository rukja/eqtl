QTL analysis pipeline
Author: Rithwik Narendra
Last Modified:  04/02/2026 RN

********************************** STEP 0: Prepare your BED file. **********************************

There is no one way of preparing a BED file from the pseudobulked data. I have included a sample R markdown which shows how one may go about doing this. The reason why there isn't a one size fits all approach is because pseudobulk data may vary based on sample names, gene names, and chromosome assembly. Nevertheless, I have included some tips to make the process easy.

First, load in the normalized NPC gene expression. Ensure that this expression has been inverse normalized. A snippet of this code is present in new_normal.R. I have pasted it below for reference

npc_df <- as.data.frame(t(as.data.frame(assay(res.proc, "NPCs"))))

int_transform <- function(x) {

  r <- rank(x, ties.method = "average")
  n <- length(r)

  u <- (r - 0.5) / n
  qnorm(u)
}

npc_int <- npc_df %>%
  as.data.frame() %>%
  mutate(across(everything(), int_transform))

Additional resources to help prepare a BED include a gene ranges file, which I have obtained from the UCSC Genome Browser for hg38 along the GENCODE V47 track. If this does not match your data, it may make sense to change the query. The code for this can be found in the markdown. The name of this file is gene_ranges.qs.

The columns of the gene ranges object must be renamed to match those in the BED file. The markdown has an example of how this was done. For appropriate indexing, ensure that the chr# format is used.

I have an additional file which is gene_locations.qs. This is not used in the current pipeline but if you feel for whatever reason that BIOMART is a better option, please feel free to use that instead. Again, reference the markdown to see how this was created and could be used.

IT IS CRITICAL that the sample names in the BED file match the sample names in the VCF in every way. There is a key that we use in the Wells Lab to rename based on old or new VCFs. QTLtools requires that sample names in the BED and VCF are in the same order. There are steps in the rest of our pipeline which perform this alignment, but it may make sense to do that upstream.

To merge the genes in the expression files with the genes in the gene ranges file, ensure that they are all in the same format (ENSG, or gene name, or ENST etc.). The file custom.txt does this for some of the NPC genes.

There are additional lines in the markdown which show how to extract condition and sample/subject ID from the columns.

The goal is for your BED file to match the BED file described in the QTLtools handbook. There are steps in the markdown which accomplish this. The key snippet is:

write.table(ctrl_bed, file = "../bigBed/ctrl.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Now add a '#' at the beginning of the header line
lines <- readLines("../bigBed/ctrl.bed")
lines[1] <- paste0("#", lines[1])  # Add '#' at the start of the first row (header)

# Write the modified lines back to the file
writeLines(lines, "../bigBed/ctrl.bed")

The desired format can be found here: https://qtltools.github.io/qtltools/ in the "Preparing input files" subheader.




********************************** STEP 1: Prepare your files and directories. **********************************

There are a number of scripts you must run for QTL analysis. Your directories must be appropriately set up for these and downstream scripts to run. This step facilitates the setup.

Create a project directory in your project/mfwells/<USERNAME> directory
1. Within that project directory, create two subdirectories called "bigVCF" and "bigBed"

2. Place the VCF file into "bigVCF" and your BED file(s) into "bigBed." Make sure your control bed is named "ctrl.bed"

3. If not already, ensure a covariates.txt file is present in your project directory. The first row should contain "ID" followed by the sample IDs. The second row should contain the covariate followed by the values:
  EX: ID  es-CHB09  es-CT4
      sex male      female

3. Move all scripts out of any directories and into the project directory. Run chmod +x *.sh to make these scripts executable. Inspect each script and ensure the notification email (if present after -M "email") is set to your email

4. cd to your project directory (if not already there) and run "qsub script_00_prepare_files.sh"

5. The structure of your files should now be (without scripts included):

[project directory]
├── beds.txt
├── covariates.txt
├── sample_names.txt
├── [bigBed]
│   ├── ctrl.bed.gz
    ├── ctrl.bed.gz.tbi
    ├── [].bed.gz.tbi
    ├── [].bed.gz
	...
├── [bigVCF]
│   ├── Cellopedia_Eggan_hESC_WGS_callset1.wellsIDs.fl3.vcf.bgz
│   ├── Cellopedia_Eggan_hESC_WGS_callset1.wellsIDs.fl3.vcf.bgz.tbi
│   ├── sample_names.txt


6. How you want to organize logs is up to you. I like to create a log folder where I can move logs into after hoffman jobs are finished.


********************************** STEP 2: Run QC and PCA on VCF. **********************************

1. Move "script_01_prepare_VCF.sh" into the bigVCF directory:
  mv script_01_prepare_VCF.sh ./bigVCF

2. Move "flashPCA.R" into the bigVCF directory:
  mv flashPCA.R ./bigVCF

3. Move "exclusion_regions_hg19.txt" into the bigVCF directory:
  mv exclusion_regions_hg19.txt ./bigVCF

4. cd into bigVCF and run "qsub ./prepare_vcf.sh [path to VCF] [MAF]"
  EX: qsub ./script_01_prepare_VCF.sh ./Cellopedia_Eggan_hESC_WGS_callset1.wellsIDs.fl3.vcf.bgz 0.20


This step performs QC on the VCF, including, in order, filtering based on BED patients and MAF frequency, HWE correction using PLINK, missingness filtering
using PLINK, filtering out individuals with low coverage USING PLINK, and then recoding back into a VCF.


More importantly, it performs PC analysis. We can inspect our PCs, which are in the output PCs.pca, by visualizing them.
This is done for you and you can inspect this in VCF_PCs.png. There are additional outputs and plots
which can facilitate PC selection.

Lastly, there are some intermediate PLINK outputs which will be used in downstream scripts.



********************************** STEP 3: Run PCA on BED. **********************************

1. cd back into the main directory
  cd ..

2. Ensure the file beds.txt is present

3. Run PCA on BED by running "qsub ./script_02_pc_bed.sh [name of bed list file] [# of genotype PCs to use] [VCF path] [whether or not to include covariates]"
  EX: qsub script_02_pc_bed.sh beds.txt 5 ./bigVCF/Cellopedia_Eggan_hESC_WGS_callset1.wellsIDs.fl3.hwe.vcf.gz N

This step performs PC and other analysis on the BED file. It, in order, calculates PCs using the script ./calculate_pcs.R, merges the PCs with the genotype PCs and covariates.txt values using ./merge_pcs.R, investigates covariate redundancy by performing correlation analysis using ./redundant_pcs.R, removes redundant pcs using ./redundant_pcs.R, reorders the pc/covariate columns so that the sample order matches that of the VCF and BED, and removes "chr" from the BED file to match the recoded VCF PLINK output.

If you would like to see which PCs were removed (if any), you can inspect the log for redundant_pcs.R in the directory: hoffman_log_Rscripts. You can additionally view the covariate correlation matrix. You can view the final covariate file in *condition*_optimal_pc_cov.txt

The covariates we choose not to use typically include sex.

********************************** STEP 4: Choose optimal # of PCs. **********************************

1. Run PC optimization using "qsub script_03_pc_optimize.sh [name of bed list file] [VCF path] [# of genotype PCs to use]"
  EX: qsub script_03_pc_optimize.sh beds.txt ./bigVCF/Cellopedia_Eggan_hESC_WGS_callset1.wellsIDs.fl3.hwe.vcf.gz 5

This step chooses the optimal number of PCs by performing CIS nominal EQTL analysis on the data using the combined covariate file created in the previous step. You can view the outputs of this analysis in the folder: pca_bed_output_*condition* by inspecting the file "pc_optimization_results.csv." The "optimal" number of PCs will be selected; the output will be in *condition*_optimal_results.txt. The optimal covariate file is now *condition*_optimal_pc_cov.txt.

Note that all VCF PCs are kept and the only PCs iterated on are the BED PCs.


********************************** STEP 5: Perform mominal and permutation analysis. **********************************

Make sure you are in your project directory. This will perform nominal and permutation analysis for all 23 chromosomes. You will do this using the script submit_all.sh

1. Create a textfile with the conditions you want to run, called "analyses.txt." On each line, write the prefix of each condition. Use "" in place of "ctrl" where "" is a blank line (do not press enter if this is the only line).
  EX:
	""
	vpa
	response

2. Run "qsub submit_all.sh"



This should create many logs, one for each chromosome, and many directories, one for each chromosome. It should create an output directory called *condition*_permutation_dir which contains the permutation results for all chromosomes.

********************************** STEP 6: Perform FDR thresholding and conditional analysis. **********************************

Make sure you are in your project directory. Make sure the script qtltools_runFDR_cis.R is present in this directory.

1. Run "qsub script_05_compute_thresholds.sh [prefix]" where prefix refers to the prefix affixed to permutation directory. For example, if my permutation directory were named ctrl_permutation_dir, my prefix would be ctrl.
  EX: qsub script_05_compute_thresholds.sh ctrl
  EX: qsub script_05_compute_thresholds.sh response


This concatenates all permutation outputs and performs FDR analysis. You can view the FDR analysis output in the hoffman_log_Rscripts log directory. It then performs conditional analysis using these thresholds.

********************************** STEP 7: Postprocessing of conditional analysis and run LD **********************************

Make sure you are in your project directory.

1. Run "qsub scruot_06_prepare_variants_plus_ld.sh [plink_prefix] [output_prefix]" where plink prefix refers to the text found before the PLINK file type identifiers. You can find this by navigating to ./bigVCF/downstream_ld_output/ (these were created by 01)
  EX: qsub scruot_06_prepare_variants_plus_ld.sh Cellopedia_Eggan_hESC_WGS_callset1.wellsIDs.fl3.hwe_plink_filtered test

This combines all conditional results, creates a list of hits, a list of re/eGenes, extracts the top conditional variants, creates a list of top/independent hits, and creates variants.txt, which is used in the next step to extract genotypes.

This script also performs within-gene LD analysis using PLINK and outputs LD matrices.


********************************** STEP 8: Extract genotypes for QTLs of interest **********************************

Make sure you are in your project directory

1. Check the variants.txt output of 06. The first column should include the variant. The second column should include the chromosome.
  EX: rs282888 chr1
      rs32922  chrX

2. Run "qsub script_06_extract_genotypes_variant_list.sh [variant_file]" where variant_file is of the form above
  EX: qsub script_06_extract_genotypes_variant_list.sh variants.txt

********************************************************************************************************************



