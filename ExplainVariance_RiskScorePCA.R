#Load Libraries
library(SummarizedExperiment)
BiocManager::install("scater")
install.packages("scater")
library(scater)
library(edgeR)
library(dplyr)
library(ggplot2)
library(matrixStats)

GeneCounts <- read.csv("pathway/GeneCount_FinalSampleSize", header=TRUE, row.names=1)
#This gene counts is not normalised here
PhenotypeInfo = read.csv("pathway/Pheno_FinalSampleSize", header=TRUE, row.names=1)
QC_Metrics_sub= read.csv("pathway/QC_Metrics_FinalSampleSize", header=TRUE, row.names=1)

REM= "SampleExcluded" #Remove participant
GeneCounts= GeneCounts[, !(colnames(GeneCounts) %in% REM)]

rownames(QC_Metrics_sub)= gsub("Structure_.*", "Structure/", rownames(QC_Metrics_sub))


rownames(PhenotypeInfo)= gsub("/", "_", rownames(PhenotypeInfo))
PhenotypeInfo[, 1] = gsub("/", "_", PhenotypeInfo[, 1])
rownames(PhenotypeInfo)= PhenotypeInfo[, 1]

GeneCounts= t(GeneCounts)
PhenotypeInfo_1 <- PhenotypeInfo[match(rownames(GeneCounts), rownames(PhenotypeInfo)), ]
QC_Metrics_sub_1 <- QC_Metrics_sub[match(rownames(PhenotypeInfo_1), rownames(QC_Metrics_sub)), ]
Phenotype_QC= cbind(PhenotypeInfo_1, QC_Metrics_sub_1)
Phenotype_QC= Phenotype_QC[!(rownames(Phenotype_QC) %in% REM), ]

GeneCounts= t(GeneCounts)

#Remove Parkinson Samples
removedSamp= Phenotype_QC$subjid[is.na(Phenotype_QC[, 6])]
print(removedSamp)
Phenotype_QC= Phenotype_QC[!is.na(Phenotype_QC[, 6]), ]
GeneCounts= GeneCounts[, !colnames(GeneCounts) %in% removedSamp]
GeneCounts=as.data.frame(GeneCounts)


all(rownames(Phenotype_QC) %in% colnames(GeneCounts))


#Prepare data structure
GeneCounts=as.matrix(GeneCounts)

#Creating summarized experiemnt
SumEx= SummarizedExperiment(
  assays=list(counts=as.matrix(GeneCounts)),
  colData= Phenotype_QC,
  rowData=data.frame(geneid= rownames(GeneCounts))
)

#Normalsie Counts
dge=DGEList(counts= GeneCounts)
colSums(dge$counts)
dge=calcNormFactors(dge)
dge$samples$norm.factors
min(dge$samples$norm.factors)
cpm_dge <- edgeR::cpm(dge)
logcpm=log2(cpm_dge+1)

#Initialise phenotyoe and qc variables 
gender=colData(SumEx)$gender.y
age=colData(SumEx)$age
ProRisk1=colData(SumEx)$post_prob_v1
PCT_PFA_Reads=colData(SumEx)$PCT_PF_READS_ALIGNED
PF_A_Bases=colData(SumEx)$PF_ALIGNED_BASES
PCT_ADAPTER=colData(SumEx)$PCT_ADAPTER
Avg_read_length=colData(SumEx)$MEAN_ALIGNED_READ_LENGTH

#Subset to look at top 500 expressed genes
meanCount=rowMeans(logcpm)
topGenes=order(meanCount, decreasing=TRUE)[1:500]
logCount_sub=logcpm[topGenes, ]

#Set up variables for the following forloop on explained variance by each variable

#Demographic/phenotype related variables
var_ex_gen= numeric(nrow(logCount_sub))
var_ex_age= numeric(nrow(logCount_sub))
var_ex_both_gen_age= numeric(nrow(logCount_sub))
var_ex_ALL= numeric(nrow(logCount_sub))
var_ex_R1= numeric(nrow(logCount_sub))


#Picard QC Metrics variables
var_ex_PCT_PFA_Reads= numeric(nrow(logCount_sub))
var_ex_PF_A_Bases= numeric(nrow(logCount_sub))
var_ex_PCT_ADAPTER= numeric(nrow(logCount_sub))
var_ex_Avg_read_length= numeric(nrow(logCount_sub))
var_ex_Avg_qc= numeric(nrow(logCount_sub))


#Calculating explained variance
for (i in 1:nrow(logCount_sub)) {
  #Gender
  lmMod_gen=lm(logCount_sub[i, ] ~ gender)
  var_ex_gen[i]=summary(lmMod_gen)$r.squared
  #Age
  lmMod_age=lm(logCount_sub[i, ] ~ age)
  var_ex_age[i]=summary(lmMod_age)$r.squared
  #All: gender, age, riskV1
  lmMod_ALL=lm(logCount_sub[i, ] ~ gender + age + ProRisk1)
  var_ex_ALL[i]=summary(lmMod_ALL)$r.squared
  #Prod Risks
  lmMod_ProRisk1=lm(logCount_sub[i, ] ~ ProRisk1)
  var_ex_R1[i]=summary(lmMod_ProRisk1)$r.squared
  #PCT_PFA_Reads
  lmMod_PCT_PFA_Reads=lm(logCount_sub[i, ] ~ PCT_PFA_Reads)
  var_ex_PCT_PFA_Reads[i]=summary(lmMod_PCT_PFA_Reads)$r.squared
  #PF_A_Bases
  lmMod_PF_A_Bases=lm(logCount_sub[i, ] ~ PF_A_Bases)
  var_ex_PF_A_Bases[i]=summary(lmMod_PF_A_Bases)$r.squared
  #PCT_ADAPTER
  lmMod_PCT_ADAPTER=lm(logCount_sub[i, ] ~ PCT_ADAPTER)
  var_ex_PCT_ADAPTER[i]=summary(lmMod_PCT_ADAPTER)$r.squared
  #Avg_read_length
  lmMod_Avg_read_length=lm(logCount_sub[i, ] ~ Avg_read_length)
  var_ex_Avg_read_length[i]=summary(lmMod_Avg_read_length)$r.squared
  #All picard qc metrics
  lmMod_All_qc=lm(logCount_sub[i, ] ~ Avg_read_length+ PCT_ADAPTER + PF_A_Bases + PCT_PFA_Reads)
  var_ex_Avg_qc[i]=summary(lmMod_All_qc)$r.squared
}


#Creating combined dataframe
varFacgen= data.frame(VarianceExplained=var_ex_gen, Variable= "Gender:__%")
GenMean= mean(varFacgen[, 1], na.rm=TRUE) * 100

varFacage= data.frame(VarianceExplained=var_ex_age, Variable= "Age:__%")
AgeMean= mean(varFacage[, 1], na.rm=TRUE) * 100

varFacboth_ALL= data.frame(VarianceExplained=var_ex_ALL, Variable= "Gender, Age & Prodromal_Risk:__%")
All_three= mean(varFacboth_ALL[, 1], na.rm=TRUE) * 100

varFac_var_ex_R1= data.frame(VarianceExplained=var_ex_R1, Variable= "Prodromal_Risk:___%")
R1Mean= mean(varFac_var_ex_R1[, 1], na.rm=TRUE) * 100

varFacPCT_PFA_Reads= data.frame(VarianceExplained=var_ex_PCT_PFA_Reads, Variable= "PCT_PFA_Reads:___%")
PCT_PFA_ReadsMean= mean(varFacPCT_PFA_Reads[, 1], na.rm=TRUE) * 100

varFacPF_A_Bases= data.frame(VarianceExplained=var_ex_PF_A_Bases, Variable= "PF_A_Bases:___%")
PF_A_BasesMean= mean(varFacPF_A_Bases[, 1], na.rm=TRUE) * 100

varFacPCT_ADAPTER= data.frame(VarianceExplained=var_ex_PCT_ADAPTER, Variable= "PCT_ADAPTER:___%")
PCT_ADAPTERMean= mean(varFacPCT_ADAPTER[, 1], na.rm=TRUE) * 100

varFacAvg_read_length= data.frame(VarianceExplained=var_ex_Avg_read_length, Variable= "Avg_read_length:___%")
Avg_read_lengthMean= mean(varFacAvg_read_length[, 1], na.rm=TRUE) * 100

varFac_ALL_QC_Metrics= data.frame(VarianceExplained=var_ex_Avg_qc, Variable= "Quality_Metrics:___%")
ALL_QC_MetricsMean= mean(varFac_ALL_QC_Metrics[, 1], na.rm=TRUE) * 100


combDF_demo= rbind(varFacgen, varFacage, varFac_var_ex_R1, varFacboth_gen_age, varFacboth_ALL)
combDF_demo_poster= rbind(varFacgen, varFacage, varFac_var_ex_R1, varFacboth_ALL, varFac_ALL_QC_Metrics)
combDF_demo_poster$Variable= factor(combDF_demo_poster$Variable, levels= c("Quality_Metrics:__%", "Gender, Age & Prodromal_Risk:__%", "Prodromal_Risk:__%", "Age:__%", "Gender:__%"))
combDF_qc= rbind(varFacPCT_PFA_Reads, varFacPF_A_Bases, varFacPCT_ADAPTER, varFacAvg_read_length, varFac_ALL_QC_Metrics)
combDF_qc$Variable= factor(combDF_qc$Variable, levels= c("ALL QC Metrics:__%", "PF_A_Bases:__%", "Avg_read_length:__%", "PCT_PFA_Reads:__%", "PCT_ADAPTER:__%"))


#Plot variance explained 
ggplot(combDF_demo_poster, aes(x = VarianceExplained*100, color=Variable)) +
  geom_density(alpha = 0.5) +
  labs(title = "Variance Explained by Clinical & Quality Metric Covariates", x = "Variance Explained (%)", y = "Density") +
  scale_fill_manual(
    values=c("blue", "red", "green", "purple", "orange"), 
    labels=c("Quality_Metrics:10.49%","Gender, Age & Prodromal_Risk:__%", "Prodromal_Risk:__%", "Age:__%", "Gender:__%"))+  theme_minimal() +
  theme(legend.title= element_blank())


ggplot(combDF_qc, aes(x = VarianceExplained*100, color=Variable)) +
  geom_density(alpha = 0.5) +
  labs(title = "Variance Explained by Picard Quality Contrl Metrics", x = "Variance Explained (%)", y = "Density") +
  scale_fill_manual(values=c("blue", "red"), labels=c("ALL QC Metrics: __%", "PF_A_Bases: __%", "Avg_read_length: __%", "PCT_PFA_Reads: __%", "PCT_ADAPTER: __%"))  +
  theme_minimal() +
  theme(legend.title= element_blank())

#To see individual genes in boxplot 
ggplot(combDF, aes(x=Variable, y= VarianceExplained, color= Variable)) +
  geom_boxplot() +
  labs(title= "Variance Explained",
       x = "Factor",
       y= "Variance Explained (R-squared)") +
  ylim(0, 0.3) +
  theme_minimal() +
  theme(legend.title=element_blank())



#Saving subsetted data for Differential Expression Analysis
PhenotypeInfo_1$post_prob_v1[is.na(PhenotypeInfo_1$post_prob_v1) & PhenotypeInfo_1$study_arm == "PD"] = 1

PhenotypeInfo_1_cl= subset(PhenotypeInfo_1, !is.na(post_prob_v1))

PhenotypeInfo_1_cl$post_prob_v2[is.na(PhenotypeInfo_1_cl$post_prob_v2) & PhenotypeInfo_1_cl$study_arm == "PD"] = 1

GeneCounts_cl= GeneCounts[, colnames(GeneCounts) %in% rownames(PhenotypeInfo_1_cl)]

GeneCounts_cl= as.data.frame(GeneCounts_cl)

write.csv(PhenotypeInfo_1_cl, "pathway/PhenotypeInfo_1_cl.csv", row.names= TRUE) 
#Subset phenotype data to only include sample ids that have a prodromal risk V1 score

write.csv(GeneCounts_cl, "pathway/GeneCounts_cl.csv", row.names= TRUE) 
#Subset gene data to only include sample ids that have a prodromal risk V1 score




#PCA: Top 500 expressed genes and post_prob V1 risk
GeneCounts= as.matrix(GeneCounts)
meanCount= rowMeans(GeneCounts)
topGenes=order(meanCount, decreasing=TRUE)[1:500]
GeneCounts_500= GeneCounts[topGenes, ]


GeneCounts_500_log=log2(GeneCounts_500+1)
pca_result= prcomp(t(GeneCounts_500_log), scale=FALSE)
pca_data= data.frame(pca_result$x)
pca_data$Risk_Score= Phenotype_QC$post_prob_v1
summary(pca_result)

ggplot(pca_data, aes(x=PC1, y=PC2, color= Risk_Score))+
  geom_point()+
  theme_minimal()+
  labs(title= "PCA on top 500 expressed genes in the RBD Cohort")

pca_subset= pca_result$x[, 1:5]
pca_df= data.frame(pca_subset)
pca_df$Risk_Score= Phenotype_QC$post_prob_v1

shapiro.test(pca_df$Risk_Score)#Not normally distributed so using kendall rank for correlation

corTest1= cor.test(pca_df$PC1, pca_df$Risk_Score, method="kendall")
corTest2= cor.test(pca_df$PC2, pca_df$Risk_Score, method="kendall")
corTest3= cor.test(pca_df$PC3, pca_df$Risk_Score, method="kendall")
corTest4= cor.test(pca_df$PC4, pca_df$Risk_Score, method="kendall")
corTest5= cor.test(pca_df$PC5, pca_df$Risk_Score, method="kendall")

p_val= c(corTest1$p.value, corTest2$p.value, corTest3$p.value, corTest4$p.value, corTest5$p.value)
print(p_val)

adjusted_p= p.adjust(p_val, method="BH")
print(adjusted_p)

pca_long= pca_df %>%
  gather(key="PC", value="PCA_value", PC1:PC5)

ggplot(pca_long, aes(x=PCA_value, y=Risk_Score))+
  geom_point(aes(color=Risk_Score), size=3)+
  geom_smooth(method="lm", color="red")+
  labs(title="Correlations between Risk Score and Principal Components of Gene Expression", x= "Principal Component Values", y="Prodromal Risk Score")+
  facet_wrap(~PC)+
  theme_minimal()+
  theme(legend.title= element_blank())









