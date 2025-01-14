#MultiQC Quality Control Check To detect Poor Quality Samples and Compare Positive Biological Samples against Negative Controls

#Load libraries
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)

general_stats= read_tsv("general_stats")

summary(general_stats)

general_stats= general_stats %>%
  mutate(SampleType= if_else(str_starts(Sample, "StartOfName_"), "Positive", "Negtaive"))

#Filter samples that are of poor quality in addition to the negative controls. Check if negative controls reflect as true controls.
poor_quality_samples= general_stats %>%
  filter(
    total_sequences < 1000000|
      percent_gc < 40 | percent_gc > 60|
      percent_duplicates >50 |
      percent_fails > 20
  )

print(poor_quality_samples)

ggplot(general_stats, aes(x= SampleType, y=total_sequences, fill= SampleType))+
  geom_boxplot()+
  labs(title= "Total Sequences", x= "Sample Type", y= "total_sequences")+
  scale_fill_manual(values = c("Positive"="skyblue", "Negative"="purple"))+
  ylim(584220, 176400114)+
  theme_minimal()+
  theme(
    legend.position="none",
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    plot.title=element_text(size=16, face="bold"),
  )


adapter_content= read_tsv("adapter_content") 

overrep_seq= read_tsv("overrep_seq")
overrep_seq= overrep_seq %>%
  mutate(SampleType= if_else(str_starts(Sample, "StartOfName_"), "Positive", "Negtaive"))

#Filter samples that are of poor quality in addition to the negative controls. Check if negative controls reflect as true controls.
poor_quality_samples= overrep_seq %>%
  filter(
    `Sum of remaining overrepresented sequences` > 10|
      `Top overrepresented sequence`>10 )

print(poor_quality_samples)

ggplot(overrep_seq, aes(x= SampleType, y=`Top overrepresented sequence`, fill= SampleType))+
  geom_boxplot()+
  labs(title= "Top overrepresented sequence", x= "Sample Type", y= "Percent")+
  scale_fill_manual(values = c("Positive"="skyblue", "Negative"="purple"))+
  ylim(0, 70)+
  theme_minimal()+
  theme(
    legend.position="none",
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    plot.title=element_text(size=16, face="bold"),
  )

per_base_n= read_tsv("per_base_n_file")

per_base_seq= read_tsv("per_base_seq_file")

per_base_seq= read_tsv("per_base_seq_file")

per_seq_gcP= read_tsv("per_seq_gcP_file")

per_seq_qual= read_tsv("per_seq_qual_file")

per_seq_counts= read_tsv("per_seq_counts_file")

Mereged_per_seq_counts$total_sequences

Mereged_per_seq_counts= merge(per_seq_counts, general_stats[, c("Sample", "total_sequences")], by="Sample", all.x = TRUE)

Mereged_per_seq_counts$DupPer=(Mereged_per_seq_counts$`Duplicate Reads`/Mereged_per_seq_counts$total_sequences)
Mereged_per_seq_counts$UniPer=(Mereged_per_seq_counts$`Unique Reads`/Mereged_per_seq_counts$total_sequences)

#Filter samples that are of poor quality in addition to the negative controls. Check if negative controls reflect as true controls.
poor_quality_samples= Mereged_per_seq_counts %>%
  filter(
    DupPer > 0.50|
      UniPer < 0.50 )

print(poor_quality_samples)

dup_levels= read_tsv("dup_levels_file")

dup_levels$Total= rowSums(dup_levels[, 11:ncol(dup_levels)], na.rm=TRUE)

dup_levels= dup_levels %>%
  mutate(SampleType= if_else(str_starts(Sample, "StartOfName_"), "Positive", "Negtaive"))

#Filter samples that are of poor quality in addition to the negative controls. Check if negative controls reflect as true controls.
poor_quality_samples= dup_levels %>%
  filter(
    Total > 25|
      `>100`>15|
      `>10k+`>15| 
      `>10`>15|
      `>500`>15|
      `>1k`>15|
      `>5k`>15|
      `>50`>15
  )

ggplot(dup_levels, aes(x= SampleType, y=Total, fill= SampleType))+
  geom_boxplot()+
  labs(title= "Duplicated Sequences", x= "Sample Type", y= "Percent")+
  scale_fill_manual(values = c("Positive"="skyblue", "Negative"="purple"))+
  ylim(0, 100)+
  theme_minimal()+
  theme(
    legend.position="none",
    axis.text=element_text(size=12),
    axis.title=element_text(size=14),
    plot.title=element_text(size=16, face="bold"),
  )

