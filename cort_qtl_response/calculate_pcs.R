#!/usr/bin/env Rscript

# R script to perform PCA and save all principal components
library(data.table)
library(tidyverse)

if (!requireNamespace("findPC", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
  }
  devtools::install_github("haotian-zhuang/findPC")
}

library(findPC)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 2) {
  stop("Exactly 2 arguments needed: <bed_file_path> <output_prefix>", call.=FALSE)
}

bed_file <- args[1]      # Path to the compressed bed file
output_prefix <- args[2] # Output file prefix for PCA results

read_bed <- function(bed_file) {
  # If your bed file is compressed
  if (grepl("\\.gz$", bed_file)) {
    cmd <- paste("zcat", bed_file)
    bed_data <- fread(cmd = cmd)
  } else {
    bed_data <- fread(bed_file)
  }
  return(bed_data)
}

# Print confirmation
cat("Reading bed file:", bed_file, "\n")
cat("Output will be saved to:", paste0(output_prefix, ".pca"), "\n")

# Read data
bed_data <- read_bed(bed_file)

# Phenotype data in this BED file begins at column 7
phenotype_data <- as.matrix(bed_data[, 7:ncol(bed_data), with=FALSE])
rownames(phenotype_data) <- bed_data$gid # The header should be read in, it starts with a #

# Center and scale the data
phenotype_data_scaled <- scale(t(phenotype_data), center=TRUE, scale=TRUE)

# Perform PCA to calculate ALL pcs
pca_result <- prcomp(phenotype_data_scaled, center=FALSE, scale=FALSE)


# Make elbow plot
sdev <- pca_result$sdev[1:ncol(phenotype_data_scaled)]

pdf("genexPC_plot.pdf")

sdev <- sort(sdev, decreasing = TRUE)
findPC(sdev = sdev, figure = TRUE)
dev.off()


p <- findPC(sdev = sdev)

# Select optimal # of PCs


n_pcs <- as.numeric(p)

# Get sample IDs (column names from your original data)
sample_ids <- colnames(phenotype_data)

# Create PCA output in QTLtools format
pca_output <- as.data.frame(pca_result$x[, 1:n_pcs, drop = FALSE])
colnames(pca_output) <- paste0("output_center_scale_PC", 1:ncol(pca_output))
pca_output <- cbind(sample_ids, pca_output)

pca_output <- t(pca_output)

# Save all PCs in QTLtools format
write.table(pca_output, file=paste0(output_prefix, ".pca"),
            quote=FALSE, sep="\t", row.names=TRUE, col.names = TRUE)

lines <- readLines(paste0(output_prefix, ".pca"))
lines[1] <- paste("SampleID", lines[1], sep = "\t")
writeLines(lines, paste0(output_prefix, ".pca"))



#Calculate variance metrics
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
cum_var <- cumsum(var_explained)

#Create a summary table
pca_summary <- data.frame(
  PC = 1:length(pca_result$sdev),
  Standard_Deviation = pca_result$sdev,
  Proportion_of_Variance = var_explained,
  Cumulative_Proportion = cum_var
)

#Save it as a TSV
write.table(pca_summary, file=paste0(output_prefix, ".pca_summary.tsv"),
            quote=FALSE, sep="\t", row.names=FALSE)

