args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
  stop("Usage: Rscript merge_covariates.R file1.txt file2.txt file3.txt inclusion output.txt")
}

# Inputs
file1 <- args[1]
file2 <- args[2]
file3 <- args[3]
inclusion_df2 <- args[4]
outfile <- args[5]

# Read files
df1 <- read.table(file1, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
df2 <- read.table(file2, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
df3 <- read.table(file3, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Merge based on common samples (column names)
common_samples <- Reduce(intersect, list(colnames(df1)[-1], colnames(df2)[-1], colnames(df3)[-1]))

if (length(common_samples) == 0) {
  stop("No common samples found across all files!")
}

# Reorder columns to match df1
reorder_columns <- function(df, common_samples) {
  c(colnames(df)[1], common_samples)
}

df1_ordered <- df1[, reorder_columns(df1, common_samples), drop = FALSE]
df2_ordered <- df2[, reorder_columns(df2, common_samples), drop = FALSE]
df3_ordered <- df3[, reorder_columns(df3, common_samples), drop = FALSE]

# Merge rows by stacking them (keeping header consistent)

if (inclusion_df2 == "Y") {
  merged_df <- rbind(df1_ordered, df2_ordered, df3_ordered[-1,])
} else {
  merged_df <- rbind(df1_ordered, df3_ordered[-1,])
}



rownames(merged_df) <- merged_df[[1]]

# Write output
write.table(merged_df, file = outfile, sep = "\t", quote = FALSE, row.names = FALSE)
