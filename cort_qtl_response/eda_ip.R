args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 1)  stop("Incorrect number of arguments")

opt_input  <- args[1] #opt_input is "./conditionals_full.txt"
temp <- read.table(opt_input,nrows = 1)

header <- FALSE
if (temp$V1 == "phe_id" || temp$V1 == "grp_id"){
  header  <- TRUE
  bwd_pval <- which(temp == "bwd_pval")
  if(length(bwd_pval) != 1)  stop("Header present but cannot find column bwd_pval")
  bwd_r_squared <- which(temp == "bwd_r_squared")
  if(length(bwd_r_squared) != 1)  stop("Header present but cannot find column bwd_r_squared")
  bwd_slope <- which(temp == "bwd_slope")
  if(length(bwd_slope) != 1)  stop("Header present but cannot find column bwd_slope")
  bwd_best_hit <- which(temp == "bwd_best_hit")
  if(length(bwd_best_hit) != 1)  stop("Header present but cannot find column bwd_best_hit")
  rank <- which(temp == "rank")
  if(length(rank) != 1)  stop("Header present but cannot find column rank")
}
D = read.table(opt_input,hea=header, stringsAsFactors=FALSE)

## Create order columns

D_ranks <- D %>%
  group_by(phe_id) %>%
  mutate(
    magnitude_bwd_slope = abs(bwd_slope),
    
    slope_rank = dense_rank(-magnitude_bwd_slope),
    rsq_rank   = dense_rank(-bwd_r_squared),
    adj_pval_rank = dense_rank(bwd_pval),
    
    score = magnitude_bwd_slope * -log10(pmax(bwd_pval, 1e-300)),
    rsq_score = magnitude_bwd_slope * -log10(pmax(bwd_pval, 1e-300)) * bwd_r_squared,
    
    score_rank = dense_rank(-score),
    rsq_score_rank = dense_rank(-rsq_score)
  )
  

ggplot(D_ranks, aes(x = as.factor(bwd_best_hit), y = score_rank)) +
  geom_boxplot()

ggplot(D_ranks, aes(x = as.factor(bwd_best_hit), y = slope_rank)) +
  geom_boxplot()

ggplot(D_ranks, aes(x = as.factor(bwd_best_hit), y = score)) +
  geom_boxplot()

## Annotation

# Reformat

variants_to_annotate <- D %>%
  dplyr::select(var_id, var_chr, var_from, var_to, phe_id, fwd_pval, fwd_slope, phe_strd) %>%
  distinct(var_id, .keep_all = TRUE) %>%
  mutate(chr_pos = paste0("chr", var_chr, ":", var_from)) %>%
  arrange(fwd_pval)

write.table(variants_to_annotate[, c("var_id", "chr_pos", "var_from", "var_to")],
            "variants_for_annotation.txt", 
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")


# Perform annotation

library(VariantAnnotation)
library(biomaRt)

#ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Create GRanges from variants
variant_gr <- GRanges(
  seqnames = paste0("chr", variants_to_annotate$var_chr),
  ranges = IRanges(start = variants_to_annotate$var_from, 
                   end = variants_to_annotate$var_to),
  strand = variants_to_annotate$phe_strd, 
  variant_id = variants_to_annotate$var_id,
  gene = variants_to_annotate$phe_id
)

BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
g.overlap_features <- IRanges:findOverlaps(variant_gr, genes(txdb), single.strand.genes.only=FALSE)
