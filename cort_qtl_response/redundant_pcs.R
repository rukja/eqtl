#!/usr/bin/env Rscript

# R script to perform correlation analysis on PCs
library(data.table)
library(tidyverse)
library(dplyr)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments is provided
if (length(args) != 2) {
  stop("Exactly 2 arguments needed: <covariate file path>, <bed base name>", call.=FALSE)
}

cov_file <- args[1]      # Path to the covariate file
bed <- args[2]

cor_threshold <- 0.9


# Print confirmation
cat("Reading covariate file:", cov_file, "\n")


# Read data
cov <- read.table(cov_file, header=TRUE, sep="", check.names=FALSE, stringsAsFactors=TRUE, row.names = 1)

cov <- as.data.frame(t(cov))

# Convert factors to numeric (and characters, if any)
cov[] <- lapply(cov, function(x) {
  if (is.factor(x)) {
    as.numeric(as.factor(x))  # or as.numeric(x) if you want factor levels as numbers
  } else if (is.character(x)) {
    as.numeric(as.factor(x))  # Converts characters to factors, then to numbers
  } else {
    x  # Leave numeric and other classes unchanged
  }
})



cov_numeric <- cov[, sapply(cov, is.numeric)]
cor_mat <- cor(cov_numeric, use = "pairwise.complete.obs")


cor_file <- paste0(bed, "_pc_cov_correlation_matrix.txt")

cat(cor_file, "\n")

write.table(cor_mat, file = cor_file, sep = "\t", row.names = TRUE, quote = FALSE)

cat("Investigating redundancy", "\n")

redundant_pairs <- which(abs(cor_mat) > cor_threshold & abs(cor_mat) < 1, arr.ind = TRUE)
redundant_df <- data.frame(
  Var1 = rownames(cor_mat)[redundant_pairs[,1]],
  Var2 = colnames(cor_mat)[redundant_pairs[,2]],
  Cor = cor_mat[redundant_pairs]
)

if (nrow(redundant_df) == 0) {
  message("No covariates found with correlation above threshold.")
} else {
  print(redundant_df)
  vars_to_remove <- unique(redundant_df$Var2[redundant_df$Var1 < redundant_df$Var2])
  message("Variables to remove (correlation > ", cor_threshold, "):")
  print(vars_to_remove)


  cov_filtered <- as.data.frame(t(cov[, !colnames(cov) %in% vars_to_remove]))
 # cov_filtered_file <- paste0(bed, "_filtered.txt")
	
  write.table(cov_filtered, file = cov_file, sep = "\t", row.names = TRUE, quote = FALSE)

}
