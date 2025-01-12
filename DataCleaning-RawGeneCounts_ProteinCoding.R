#Initialise packages
install.packages("BiocManager")
BiocManager::install("rtracklayer")
library(rtracklayer)
BiocManager::install("GenomicRanges")
library(GenomicRanges)
library(tidyr)
library(tibble)
library(dplyr)
library(tibble)

#Load Data
Raw_GeneCounts <- read.csv("pathway/Raw_GeneCounts.csv", header=TRUE, row.names=1)
Prodromal_Data = read.csv("pathway/Prodromal_Data.csv", header=TRUE, row.names=1)
Quality_Data = read.csv("pathway/Quality_Data.csv", header=TRUE, row.names=1)
Prodromal_Quality_Data= cbind(Prodromal_Data, Quality_Data)

#remove poor quality samples and make sure Gene Count Samples are matched to Prodromal_Quality_Data Samples 
REM= c("SampleNames", "....") 
Raw_GeneCounts= Raw_GeneCounts[, !(colnames(Raw_GeneCounts) %in% REM)]
rownames(Prodromal_Quality_Data)=Prodromal_Quality_Data$SampleName
Prodromal_Quality_Data= Prodromal_Quality_Data[match(colnames(Raw_GeneCounts), rownames(Prodromal_Quality_Data)), ]
Raw_GeneCounts=as.data.frame(Raw_GeneCounts)

#Access gene names from Reference Annotation file to subset the Raw_GeneCounts to only contain Protein Coding Genes
gtf_1= "pathway/reference_annotation"
gtf_1dat= import(gtf_1)
gtf_1df= as.data.frame(gtf_1dat)
gtf_1dat_sub <- gtf_1df %>% select(seqnames, Geneid, gene_name, gene_type)

ProteinCoding_dat= gtf_1dat_sub %>% #Subset the reference file to only contain protein coding genes
  filter(gene_type == "protein_coding")

ProteinCoding_dat <- ProteinCoding_dat %>% select(Geneid, gene_name)

ProteinCoding_dat= ProteinCoding_dat%>%
  distinct(Geneid, .keep_all = TRUE)

MultNames= ProteinCoding_dat %>%
  group_by(Geneid) %>%
  summarize(unique_names = n_distinct(gene_name)) %>%
  filter(unique_names > 1)

print(MultNames)

#Subset Raw_GeneCounts to include only protein coding genes 
Raw_GeneCounts <- Raw_GeneCounts[rownames(Raw_GeneCounts) %in% ProteinCoding_dat$Geneid, ]
Raw_GeneCounts$Geneid=rownames(Raw_GeneCounts)
Raw_GeneCounts <- Raw_GeneCounts %>%
  select(Geneid, everything())

Raw_GeneCounts_Symbols = merge(Raw_GeneCounts, ProteinCoding_dat, by="Geneid") Merege the associated gene symbol with the Geneid
Raw_GeneCounts_Symbols <- Raw_GeneCounts_Symbols %>%
  select(gene_name, everything())
Raw_GeneCounts_Symbols= data.frame(Raw_GeneCounts_Symbols)

# Calculate average gene count across samples
Raw_GeneCounts_Symbols <- Raw_GeneCounts_Symbols %>%
  # Select only the sample columns
  mutate(across(starts_with("__"), as.numeric)) %>% # The "" is where the beginning term of sample name is writtern
  rowwise() %>%
  mutate(avgCount = mean(c_across(starts_with("__")), na.rm = TRUE)) %>% # The "" is where the beginning term of sample name is writtern
  ungroup()

Raw_GeneCounts_Symbols <- Raw_GeneCounts_Symbols %>%
  select(gene_name, Geneid, avgCount, everything())

# Filter gene names with more than one geneid
RepDat <- Raw_GeneCounts_Symbols %>%
  group_by(gene_name) %>%
  filter(n() > 1) %>%  # Retain gene_names with more than one gene_id
  arrange(gene_name, Geneid) %>%  # To list gene_ids consecutively within each gene name
  ungroup()  

RepDat_sub <- RepDat %>%
  arrange(gene_name, Geneid) %>%  # Arrange by genename and geneid
  select(gene_name, Geneid, avgCount)

#Identify the gene_id with lower expression
compareID <- RepDat_sub %>%
  group_by(gene_name) %>%
  arrange(gene_name, desc(avgCount)) %>%
  mutate(low_geneid = ifelse(avgCount == min(avgCount), Geneid, NA)) %>%
  ungroup()

compareID <- compareID %>%
  mutate(zero_avg = ifelse(avgCount == 0, "Yes", "No"))

comparison_table <- compareID %>%
  select(gene_name, Geneid, avgCount, low_geneid, zero_avg) %>%
  distinct()

low_geneid_list <- comparison_table %>%
  filter(!is.na(low_geneid)) %>%  # Exclude rows where high_geneid is NA
  pull(low_geneid) 

compareID_zeros= compareID %>%
  group_by(gene_name) %>%
  filter(all(!is.na(low_geneid))) %>%
  ungroup()

compareID_zeros <- compareID_zeros %>%
  group_by(gene_name) %>%
  arrange(gene_name, desc(avgCount)) %>%
  mutate(low_geneid = ifelse(avgCount == min(avgCount),ifelse(row_number()==1, Geneid, NA), Geneid)) %>%
  ungroup()

compareID_zeros <- compareID_zeros %>%
  select(gene_name, Geneid, avgCount, low_geneid, zero_avg) %>%
  distinct()

low_geneid_list_zero <- compareID_zeros %>%
  filter(!is.na(low_geneid)) %>%  # Exclude rows where high_geneid is NA
  pull(low_geneid) 

Duplicates_genename=low_geneid_list[!low_geneid_list %in% low_geneid_list_zero] #Defines the list of genes to be removed if they are matched to the same gene symbol as another gene id, however, these gene ids have lower average expression so they will be rmeoved

#Removing Duplicates to obatin a cleaned dataset where all gene names are present once
Raw_GeneCounts_Symbols= as.data.frame(Raw_GeneCounts_Symbols)
rownames(Raw_GeneCounts_Symbols)= Raw_GeneCounts_Symbols$Geneid
Raw_GeneCounts_Symbols=Raw_GeneCounts_Symbols[!rownames(Raw_GeneCounts_Symbols) %in% Duplicates_genename, ]
rownames(Raw_GeneCounts_Symbols)= NULL
rownames(Raw_GeneCounts_Symbols)= Raw_GeneCounts_Symbols$gene_name

Raw_GeneCounts_Symbols_Cleaned=Raw_GeneCounts_Symbols[ , -c(1, 2, 3)] #Save file for Differential Gene Expression Analysis
