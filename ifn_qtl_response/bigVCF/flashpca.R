#!/usr/bin/env Rscript
# Load necessary libraries
library(data.table)
library(dplyr)
library(flashpcaR)
library(GGally)
library(findPC)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Please provide the PLINK prefix as argument.\nUsage: Rscript <plink_prefix>")
}

fn <- args[1]          # first argument = your PLINK file prefix
# Check if files exist
if (!file.exists(paste0(fn, ".bed")) || 
    !file.exists(paste0(fn, ".bim")) || 
    !file.exists(paste0(fn, ".fam"))) {
  stop("Required PLINK files not found!")
}

# Get number of samples
fam <- read.table(paste0(fn, ".fam"))
n_samples <- nrow(fam)

n_pcs <- min(15, n_samples - 1)

# Run flashpca with try to catch errors
f2 <- try(flashpca(fn, ndim = n_pcs))

# Check if flashpca succeeded
if (inherits(f2, "try-error")) {
  cat("Error running flashpca. Trying with minimal parameters...\n")
  # Try with minimal settings
  f2 <- flashpca(fn, ndim = 2)  # Try just 2 PCs
}

# Process results
if (!inherits(f2, "try-error") && is.list(f2) && !is.null(f2$vectors)) {
  
  
  
  pcs <- as.data.frame(f2$vectors)
  # Only keep columns up to the number of PCs actually computed
  pcs <- pcs[, 1:ncol(f2$vectors), drop = FALSE]
  
  rownames(pcs) <- read.table("sample_names.txt")$V1

  pcs <- t(pcs)
  
  rownames(pcs) <- paste0("PC", seq_len(nrow(pcs)))
  
  # Write output
  write.table(pcs, "PCs.pca", quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
  
  lines <- readLines("PCs.pca")
  lines[1] <- paste("SampleID", lines[1], sep = "\t")
  writeLines(lines, "PCs.pca")
  
  cat("PCA completed successfully. Output written to PCs.pca\n")

  if (!is.null(f2$values)) {
    # Calculate variance explained
    var_explained <- f2$values^2 / sum(f2$values^2) * 100
    

    sdev <- sqrt(f2$values)
    pdf("VCF_genexPC_plot.pdf")
    sdev_sorted <- sort(sdev, decreasing = TRUE)
    max_pcs <- length(sdev_sorted)
    search_range <- min(max_pcs, 50)
    findPC(sdev = sdev_sorted, number = search_range, figure = TRUE)
    dev.off()

    print(search_range)
    
    p <- findPC(sdev = sdev_sorted, number = search_range)
    n_pcs <- as.numeric(p)
    
    cat("\nOptimal number of PCs (determined by findPC):", n_pcs, "\n")
  } else {
    cat("Warning: Eigenvalues not available for elbow plot\n")
  }


} else {
  stop("Failed to compute PCs. Check your input data.")
}

a <- read.table("PCs.pca", header = TRUE, row.names = 1)
pca_t <- as.data.frame(t(as.matrix(a)))
pca_t$SampleID <- rownames(pca_t)

