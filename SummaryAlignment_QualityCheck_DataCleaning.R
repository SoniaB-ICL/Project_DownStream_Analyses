#Initialise packages
install.packages("ggplot2")
library(ggplot2)
library(tibble)
library(reshape2)
library(dplyr)
library(edgeR)
library(ggfortify)

#Initialise data
PairDS= read.csv("pathway/Pair_quality_metrics.csv", row.names =1) #Uploaded paired quality metrics data

#Removing columns that have all zeros or constants 
PairDS= PairDS[, -1]
Remove_Col= c(3, 4, 11, 15, 16, 17, 18, 19, 20, 23, 24, 25, 26, 27, 31)
PairDS= PairDS[, -Remove_Col]
PairDS=rownames_to_column(PairDS, var='SampleID')


#Fixing incorrect sample ids structures
PairDS= PairDS %>%
  mutate(SampleID = case_when(
    SampleID == "____" ~ "____b", 
    SampleID == "____" ~ "______a",
    TRUE ~ SampleID
  ))

#Detect Sample IDs Appearing Twice
PairDS$CommonSamp <- sub("([a-zA-Z0-9_]+)[ab](_S\\d+)$", "\\1", PairDS$SampleID)
RepID <- unique(PairDS$CommonSamp[duplicated(PairDS$CommonSamp) | duplicated(PairDS$CommonSamp, fromLast = TRUE)])
RepDat <- PairDS[PairDS$CommonSamp %in% RepID, ]
print(RepDat)

RepDat$PCT_ADAPTER <- as.numeric(RepDat$PCT_ADAPTER)
RepDat <- RepDat %>% arrange(CommonSamp, SampleID)


#Determine which of the two samples from the same participant should be included by comparing QC Metrics
BestRep <- RepDat %>%
  group_by(CommonSamp) %>%
  summarise(
    Sample_1 = SampleID[1],
    Sample_2 = SampleID[2],
    PCT_ADAPTER_1 = PCT_ADAPTER[1],
    PCT_ADAPTER_2 = PCT_ADAPTER[2],
    PCT_ADAPTER_Better = ifelse(PCT_ADAPTER[1] <= PCT_ADAPTER[2], SampleID[1], SampleID[2]),
    TOTAL_READS_1 = TOTAL_READS[1],
    TOTAL_READS_2 = TOTAL_READS[2],
    TOTAL_READS_Better = ifelse(TOTAL_READS[1] >= TOTAL_READS[2], SampleID[1], SampleID[2]),
    PF_READS_1 = PF_READS[1],
    PF_READS_2 = PF_READS[2],
    PF_READS_Better = ifelse(PF_READS[1] >= PF_READS[2], SampleID[1], SampleID[2]),
    PF_READS_ALIGNED_1 = PF_READS_ALIGNED[1],
    PF_READS_ALIGNED_2 = PF_READS_ALIGNED[2],
    PF_READS_ALIGNED_Better = ifelse(PF_READS_ALIGNED[1] >= PF_READS_ALIGNED[2], SampleID[1], SampleID[2]),
    PCT_PF_READS_ALIGNED_1 = PCT_PF_READS_ALIGNED[1],
    PCT_PF_READS_ALIGNED_2 = PCT_PF_READS_ALIGNED[2],
    PCT_PF_READS_ALIGNED_Better = ifelse(PCT_PF_READS_ALIGNED[1] >= PCT_PF_READS_ALIGNED[2], SampleID[1], SampleID[2]),
    PF_ALIGNED_BASES_1 = PF_ALIGNED_BASES[1],
    PF_ALIGNED_BASES_2 = PF_ALIGNED_BASES[2],
    PF_ALIGNED_BASES_Better = ifelse(PF_ALIGNED_BASES[1] >= PF_ALIGNED_BASES[2], SampleID[1], SampleID[2]),
    PF_HQ_ALIGNED_BASES_1 = PF_HQ_ALIGNED_BASES[1],
    PF_HQ_ALIGNED_BASES_2 = PF_HQ_ALIGNED_BASES[2],
    PF_HQ_ALIGNED_BASES_Better = ifelse(PF_HQ_ALIGNED_BASES[1] >= PF_HQ_ALIGNED_BASES[2], SampleID[1], SampleID[2]),
    PF_HQ_ALIGNED_Q20_BASES_1 = PF_HQ_ALIGNED_Q20_BASES[1],
    PF_HQ_ALIGNED_Q20_BASES_2 = PF_HQ_ALIGNED_Q20_BASES[2],
    PF_HQ_ALIGNED_Q20_BASES_Better = ifelse(PF_HQ_ALIGNED_Q20_BASES[1] >= PF_HQ_ALIGNED_Q20_BASES[2], SampleID[1], SampleID[2]),
    PF_MISMATCH_RATE_1 = PF_MISMATCH_RATE[1],
    PF_MISMATCH_RATE_2 = PF_MISMATCH_RATE[2],
    PF_MISMATCH_RATE_Better = ifelse(PF_MISMATCH_RATE[1] <= PF_MISMATCH_RATE[2], SampleID[1], SampleID[2]),
    PF_HQ_ERROR_RATE_1 = PF_HQ_ERROR_RATE[1],
    PF_HQ_ERROR_RATE_2 = PF_HQ_ERROR_RATE[2],
    PF_HQ_ERROR_RATE_Better = ifelse(PF_HQ_ERROR_RATE[1] <= PF_HQ_ERROR_RATE[2], SampleID[1], SampleID[2]),
    PF_INDEL_RATE_1 = PF_INDEL_RATE[1],
    PF_INDEL_RATE_2 = PF_INDEL_RATE[2],
    PF_INDEL_RATE_Better = ifelse(PF_INDEL_RATE[1] >= PF_INDEL_RATE[2], SampleID[1], SampleID[2]),
    READS_ALIGNED_IN_PAIRS_1 = READS_ALIGNED_IN_PAIRS[1],
    READS_ALIGNED_IN_PAIRS_2 = READS_ALIGNED_IN_PAIRS[2],
    READS_ALIGNED_IN_PAIRS_Better = ifelse(READS_ALIGNED_IN_PAIRS[1] >= READS_ALIGNED_IN_PAIRS[2], SampleID[1], SampleID[2]),
    PCT_CHIMERAS_1 = PCT_CHIMERAS[1],
    PCT_CHIMERAS_2 = PCT_CHIMERAS[2],
    PCT_CHIMERAS_Better = ifelse(PCT_CHIMERAS[1] <= PCT_CHIMERAS[2], SampleID[1], SampleID[2]),
    PCT_SOFTCLIP_1 = PCT_SOFTCLIP[1],
    PCT_SOFTCLIP_2 = PCT_SOFTCLIP[2],
    PCT_SOFTCLIP_Better = ifelse(PCT_SOFTCLIP[1] <= PCT_SOFTCLIP[2], SampleID[1], SampleID[2]),
    AVG_POS_3PRIME_SOFTCLIP_LENGTH_1 = AVG_POS_3PRIME_SOFTCLIP_LENGTH[1],
    AVG_POS_3PRIME_SOFTCLIP_LENGTH_2 = AVG_POS_3PRIME_SOFTCLIP_LENGTH[2],
    AVG_POS_3PRIME_SOFTCLIP_LENGTH_Better = ifelse(AVG_POS_3PRIME_SOFTCLIP_LENGTH[1] <= AVG_POS_3PRIME_SOFTCLIP_LENGTH[2], SampleID[1], SampleID[2]),
    .groups = 'drop'
  ) 

#Inspect output to determine the sample that should be removed due to poorer QC metric quality

#Clean dataset so it only includes unique sample ids
PairDS = PairDS %>%
  filter(!SampleID %in% c("Names of sample ids to remove"))

#Save Updated dataset
write.csv(PairDS, "pathway/PairDS_sub.csv", row.names= FALSE) 
PairDS_sub= read.csv("pathway/PairDS_sub.csv", row.names =1)
PairDS_sub= PairDS_sub[, -18]


#Detecting outliers based on insepct quality mertics based on outlier thresholds 
#Threshold detection
outliers = PairDS_sub[PairDS_sub$AVG_POS_3PRIME_SOFTCLIP_LENGTH > 12, ] #Keep replacing with each variable and specific threshold
print(outliers)

#Outlier Detected based on thresholds: 
#Total read- need to follow up
#PF_READS- none
#PF_READS_ALIGNED- none
#PCT_PF_READS_ALIGNED:Threshold<0.75, 
#PF_ALIGNED_BASES- need to follow up
#PF_HQ_ALIGNED_READS - none
#PF_HQ_ALIGNED_BASES - none
#PF_HQ_ALIGNED_Q20_BASES - none
#PF_MISMATCH_RATE - none
#PF_HQ_ERROR_RATE: Threshold > 0.005
#PF_INDEL_RATE - none
#Mean aligne dread length- need to follow up
#READS_ALIGNED_IN_PAIRS - none
#PCT_CHIMERAS - none
#PCT_ADAPTER - none
#PCT_SOFTCLIP - none
#AVG_POS_3PRIME_SOFTCLIP_LENGTH - none above 12


PairDS_sub_1 = PairDS_sub[!(rownames(PairDS_sub) == "Outlier_Participant"), ] #Removed outlier participant to prevent data skews and using poor quality sample reads

#Save Updated dataset
write.csv(PairDS_sub_1, "pathway/PairDS_sub_1.csv", row.names= TRUE) 


#Boxplot: Visualising Outliers
ggplot(PairDS_sub, aes(x = factor(1), y = PCT_ADAPTER)) + 
  geom_boxplot() + 
  theme(axis.title.x = element_blank())


#Further outlier detection
PairDS_sub$z_score_PCT_ADAPTER <- scale(PairDS_sub$PCT_ADAPTER)
PairDS_sub$outlier <- ifelse(abs(PairDS_sub$z_score_PCT_ADAPTER) > 3, "Outlier", "Not Outlier")


#PCA on quaity metrics
pca <- prcomp(PairDS_sub[,c("PCT_ADAPTER", "PCT_CHIMERAS", "PF_READS", "TOTAL_READS")], scale = TRUE)
pca_data <- as.data.frame(pca$x)
ggplot(pca_data, aes(x = PC1, y = PC2)) + geom_point()

#IOR
IQR_value <- IQR(PairDS_sub$TOTAL_READS)
lower_bound <- quantile(PairDS_sub$TOTAL_READS, 0.25) - 1.5 * IQR_value
upper_bound <- quantile(PairDS_sub$TOTAL_READS, 0.75) + 1.5 * IQR_value
outliers <- PairDS_sub[PairDS_sub$TOTAL_READS > upper_bound, ]


#Scale data to account for different ranges/units 
Scale_DS= scale(PairDS_sub)
Red_Scaled= Scale_DS[, apply(Scale_DS, 2, var) != 0] #Removes columns with 0 variance

#Structure data in the format of data frame for ease in following analysis
Red_Scaled=as.data.frame(Red_Scaled)
Red_Scaled=rownames_to_column(Red_Scaled, var='SampleID')
melt_RS <- melt(Red_Scaled) # Convert metadata to long format
Red_Scaled <- dcast(melt_RS, SampleID ~ variable, value.var = "value") #Convert data into a wide format

#Box plot showing distribution of datapoints for each QC metric
ggplot(melt_RS, aes(x = variable, y = value)) +
  geom_boxplot() +
  labs(title = "Distribution of QC Metrics") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

#Plotting the clustering of the samples based on QC Metrics
Red_Scaled_q= as.matrix(Red_Scaled[,-1])
distance_m=dist(Red_Scaled_q)
HClust= hclust(distance_m)
plot(HClust, labels=Red_Scaled$SampleID, main='Clustering by QC Metrics')

#PCA on Samples with QC Metrics as Factors
pcaQC <- prcomp(Red_Scaled_q, center= FALSE, scale. = FALSE)
print(pcaQC)
summary(pcaQC) #Identfied main contributions

#Visualise sample clustering using first to PCs based on QC metrics- this is done to idnetofy outliers/skews
plot(pcaQC$x[, 1], pcaQC$x[, 2], 
     xlab = "PC1", ylab = "PC2", 
     main = "PCA of Samples Based on QC Metrics")








