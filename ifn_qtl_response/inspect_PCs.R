# First, inspect VCF PCs and clustering
library(GGally)


vcf_pcs <- read.table("~/project-mfwells/qtl_one/bigVCF/PCs.pca", header = TRUE, row.names = 1, stringsAsFactors = FALSE)

vcf_pcs <- as.matrix(vcf_pcs)
storage.mode(vcf_pcs) <- "numeric"


pca_t <- as.data.frame(t(as.matrix(vcf_pcs)))
pca_t$SampleID <- rownames(pca_t)

ggpairs(pca_t, columns = grep("^PC", names(pca_t)))
