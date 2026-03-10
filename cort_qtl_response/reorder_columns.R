# Load VCF sample order
vcf_samples <- readLines("vcf_samples.txt")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript reorder_files.R <BED file> <BED output>")
}
bed_file <- args[1]
bed_output <- args[2]



bed <- read.table(bed_file, header=TRUE, sep="\t", check.names=FALSE, comment.char = "")

print(colnames(bed))

non_sample_cols <- c("#chr", "start", "end", "pid", "gid", "strand")  # Change if needed
sample_cols <- setdiff(colnames(bed), non_sample_cols)

# Reorder BED sample columns
bed_reordered <- cbind(
  bed[, non_sample_cols, drop=FALSE],
  bed[, vcf_samples, drop=FALSE]
)

# Write reordered BED

write.table(bed_reordered, file=bed_output, sep="\t", quote=FALSE, row.names=FALSE)

