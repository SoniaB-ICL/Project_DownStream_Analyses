#Cleaning phenotype data
RBD_dem = read.csv("/pathway/RBD_Data")
PD_dem = read.csv("/pathway/PD_Data")
meta = read.csv("/pathway/Data/Heinzel2019")

meta <- meta %>%
  mutate(visit_date = as.Date(visit_date, format = "%d%b%Y"))  # Convert date to Date class

meta_sub <- meta %>%
  group_by(subjid) %>%                             
  filter(visit_date == max(visit_date)) %>%               
  ungroup() 

colnames(PD_dem)[colnames(PD_dem)=='subjid'] = "subjid"
colnames(RBD_dem)[colnames(RBD_dem)=='subjid'] = "subjid"
combined_data <- bind_rows(PD_dem, RBD_dem)
combined_data$subjid= as.character(combined_data$subjid)
meta_sub$subjid=as.character(meta_sub$subjid)
merged_data <- full_join(meta_sub, combined_data, by = "subjid")

write.csv(merged_data, file = file.path(pathway, "cleaned_meta_DS.csv"), row.names = FALSE)


#Aligning Gene and Phenotype data to a subset dataset for both datasets so all samples are uniquely present in both samples and each participant only has only sample id
UniSamp = read.csv("pathway/PairDS_sub_1.csv")
rownames(UniSamp)=UniSamp$X
UniSamp=UniSamp[, -1]

#Subseting gene count data to only have samples from UniSamp
Pre_exprTable <- read.csv("pathway/RBD_gene_count_matrix.csv", header=TRUE, row.names=1)
Pre_exprTable <- subset(Pre_exprTable, select = -OuterSample)
Pre_exprTable <- subset(Pre_exprTable, select = -OuterSample)
colnames(Pre_exprTable)[colnames(Pre_exprTable) == "Sample"] <- "Sampleb"
commonIDS <- rownames(UniSamp)
matchSamp <- colnames(Pre_exprTable) %in% commonIDS
subset_GC_Dat <- Pre_exprTable[, matchSamp]

library(dplyr)

merged_data= merged_data %>%
  mutate(IDSAMPLE= gsub("/", "_", subjid))

colnames(subset_GC_Dat)= gsub("Structure_1.*", "Structure/2", colnames(subset_GC_Dat))
commonID= intersect(colnames(subset_GC_Dat), merged_data$IDSAMPLE)

meta_sub_sub= merged_data %>% filter(IDSAMPLE %in% commonID)
subset_GC_Dat_sub= subset_GC_Dat[, commonID, drop=FALSE]
subset_GC_Dat_sub=as.data.frame(subset_GC_Dat_sub)

meta_sub <- read.csv("pathway/cleaned_meta_DS.csv", header=TRUE)

onlyGC= setdiff(colnames(subset_GC_Dat), meta_sub$IDSAMPLE)
print(onlyGC)


#Saving data subsets that now contains samples only found in both phenotyp and clinical. each sample is unique and no two samples belong to same participant
write.csv(subset_GC_Dat_sub, "pathway/GeneCount_FinalSampleSize.csv", row.names= TRUE) 
GC_DATA= read.csv("pathway/GeneCount_FinalSampleSize.csv", header=TRUE, row.names=1)
write.csv(meta_sub_sub, "pathway/Pheno_FinalSampleSize.csv", row.names= TRUE) 

QC_Metrics= read.csv("pathway/PairDS_sub_1.csv", header=TRUE, row.names=1)
QC_Metrics_sub= QC_Metrics[colnames(subset_GC_Dat_sub),  , drop=FALSE]
write.csv(QC_Metrics_sub, "pathway/QC_Metrics_FinalSampleSize.csv", row.names= TRUE) 


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
meta = read.csv("pathway/Phenotype_QC.csv")
metaPD = read.csv("pathway/PD_Data")

metaPD <- metaPD %>%
  mutate(visit_date = as.Date(visit_date, format = "%d%b%Y"))  # Convert date to Date class

metaPD <- metaPD %>%
  group_by(subjid) %>%                             # Group by participant ID
  filter(visit_date == max(visit_date)) %>%                # Keep the row with the most recent date for each group
  ungroup() 

metaPD= metaPD %>%
  mutate(subjid= gsub("/", "_", subjid)) %>%
  select(subjid, age, everything())

metaPD$visit_date= as.character(metaPD$visit_date)

metaUP= meta %>%
  left_join(metaPD, by="subjid") %>%
  mutate(
    age= ifelse(is.na(age.x), age.y, age.x),
    visit= ifelse(is.na(visit.x), visit.y, visit.x),
    visit_date= ifelse(is.na(visit_date.x), visit_date.y, visit_date.x), 
    disease_duration_diag= ifelse(is.na(disease_duration_diag.x), disease_duration_diag.y, disease_duration_diag.x),
    disease_duration_onset= ifelse(is.na(disease_duration_onset.x), disease_duration_onset.y, disease_duration_onset.x)
  ) %>%
  select(subjid, visit, visit_date, age, disease_duration_diag, disease_duration_onset, everything(), - visit.x, -visit_date.x, -age.x, -disease_duration_diag.x, -disease_duration_onset.x, -visit.y, -visit_date.y, -age.y, -disease_duration_diag.y, -disease_duration_onset.y)


write.csv(metaUP, "pathway/Phenotype_QC_UP.csv", row.names= TRUE) 

