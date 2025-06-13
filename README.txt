QTL analysis pipeline
Author: Rithwik Narendra
Last Modified:  6/04/2025 RN

********************************** STEP 0: Prepare your BED file. **********************************

There is no one way of preparing a BED file from the pseudobulked data. I have included a sample R markdown which shows how one may go about doing this. The reason why there isn't a one size fits all approach is because pseudobulk data may vary based on sample names, gene names, and chromosome assembly.

First, load in the normalized gene expression. Ensure that this expression has been inverse normalized. Note that Dreamlet's process assays normalizes expression based on covariates but does not inverse normalize the data. Both have to be done. Normalizing the data here means that you DO NOT need to include --normal in ANY qtltools commands.

A BED file contains the following: chromosome, gene and phenotype name (gid, pid) which are the same for our data, start, end, strand, and the samples. All values are gene expression. Therefore, we need a file which provides us with the locations of each gene. We also need to ensure that the gene name in the locations file is in the SAME FORMAT as the gene name in our expression data. A helper file to facilitate this bridge may therefore also be needed.

Based on the format of your pseudobulk object, the subject and condition name may be concatenated in some way in the sample name. We need to extract the subject and the condition separately because the subject is what will match with the vCF (the VCF is the same across conditions).

Resources:

To help prepare for your BED, I include a gene ranges file, which I have obtained from the UCSC Genome Browser for hg38 along the GENCODE V47 track. If this does not match your data, it may make sense to change the query. The code for this can be found in the markdown. The name of this file is gene_ranges.qs.

The columns of the gene ranges object must be renamed to match those in the BED file. The markdown has an example of how this was done. For appropriate indexing, ensure that the chr# format is used.

I have an additional file which is gene_locations.qs. This is not used in the current pipeline but if you feel for whatever reason that BIOMART is a better option, please feel free to use that instead. Again, reference the markdown to see how this was created and could be used.

IT IS CRITICAL that the sample names in the BED file match the sample names in the VCF in every way. There is a key that we use in the Wells Lab to rename based on old or new VCFs. QTLtools requires that sample names in the BED and VCF are in the same order. There are steps in the rest of our pipeline which perform this alignment, but it may make sense to do that upstream depending on what your specific pipeline looks like.

To merge the genes in the expression files with the genes in the gene ranges file, ensure that they are all in the same format (ENSG, or gene name, or ENST etc.). The file custom.txt does this for some of the NPC genes.

There are additional lines in the markdown which show how to extract condition and sample/subject ID from the columns.

The goal is for your BED file to match the BED file described in the QTLtools handbook. There are steps in the markdown which accomplish this. This includes the # at the start and the delimited nature of the file. The key snippet is:

write.table(ctrl_bed, file = "../bigBed/ctrl.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

# Now add a '#' at the beginning of the header line
lines <- readLines("../bigBed/ctrl.bed")
lines[1] <- paste0("#", lines[1])  # Add '#' at the start of the first row (header)

# Write the modified lines back to the file
writeLines(lines, "../bigBed/ctrl.bed")

The desired format can be found here: https://qtltools.github.io/qtltools/ in the "Preparing input files" subheader.




