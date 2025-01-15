# Load libraries
library(data.table)
install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
library(dplyr)
library(ggplot2)
library(tidyr)

# List of subdirectories
subdirs <- list.dirs(dir_path, full.names = TRUE, recursive = FALSE)
count_list <- list()
# Loop through each sample to get featureCounts files
for (subdir in subdirs) {
  feature_file <- file.path(subdir, paste0(basename(subdir), ".txt"))
  if (file.exists(feature_file)) {
    # featureCounts reads
    data <- fread(feature_file, skip = 1) 
    # Geneid & count columns retained
    gene_counts <- data[, .(Geneid, Count = data[[7]])]
    # Rename: count column --> sample name
    sample_name <- basename(subdir)
    setnames(gene_counts, "Count", sample_name)
    count_list[[sample_name]] <- gene_counts
  }
}
# Merge data for all geneids across all samples
count_matrix <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_list)
# NA replaced to 0 for missing genes
count_matrix[is.na(count_matrix)] <- 0


write.csv(count_matrix, file = file.path(dir_path, "RBD_gene_count_matrix.csv"), row.names = FALSE)
print("Gene count matrix created and saved as gene_count_matrix.csv")

count_matrix <- read.csv("pathway/RBD_gene_count_matrix.csv")

# View the first rows 
head(count_matrix)
