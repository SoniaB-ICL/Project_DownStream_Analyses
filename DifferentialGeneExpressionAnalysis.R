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
install.packages("DT")
install.packages("data.table")
library(DT)
library(data.table)
library(limma)
library(edgeR)
install.packages("statmod")
library(statmod)
install.packages("nlme")
install.packages("quantreg")
install.packages("ppcor")
install.packages("MASS")
library(nlme)
library(quantreg)
library(ppcor)

#Load Data
Prodromal_Quality_Data= cbind(Prodromal_Data, Quality_Data)
GeneCount=Raw_GeneCounts_Symbols_Cleaned
all(colnames(GeneCount) == rownames(Prodromal_Quality_Data)) #Ensure all sample names in the GeneCounts Data are in the Prodromal_Quality_Data 

#Remove samples that are not in the RBD arm
Prodromal_Quality_Data= Prodromal_Quality_Data %>%
  filter(study_arm !="__") 

GeneCount <- GeneCount[ , match(rownames(Prodromal_Quality_Data), colnames(GeneCount))] #Ensure sample names are matched in the GeneCounts Data and Prodromal_Quality_Data 

#Ensure predictor and covariate variables are read as numerical or categorical depending on variable type
Prodromal_Quality_Data$post_prob_v1= as.numeric(Prodromal_Quality_Data$post_prob_v1)
Prodromal_Quality_Data$age= as.numeric(Prodromal_Quality_Data$age)
Prodromal_Quality_Data$gender.y= factor(Prodromal_Quality_Data$gender.y, levels = c("male", "female"))
Prodromal_Quality_Data$PF_ALIGNED_BASES= as.numeric(Prodromal_Quality_Data$PF_ALIGNED_BASES)
Prodromal_Quality_Data$PCT_PF_READS_ALIGNED= as.numeric(Prodromal_Quality_Data$PF_ALIGNED_BASES)
Prodromal_Quality_Data$PCT_ADAPTER= as.numeric(Prodromal_Quality_Data$PCT_ADAPTER)
Prodromal_Quality_Data$MEAN_ALIGNED_READ_LENGTH= as.numeric(Prodromal_Quality_Data$MEAN_ALIGNED_READ_LENGTH)

#Variable list defined
vars=c("gender.y", "age" , "PF_READS_ALIGNED" , "PF_ALIGNED_BASES" , "PCT_ADAPTER" , "MEAN_ALIGNED_READ_LENGTH" , "post_prob_v1")

#Detect if any variables are highly correlated
Prodromal_Quality_Data= Prodromal_Quality_Data[ , colnames(Prodromal_Quality_Data) %in% vars]
cor_matrix <- cor(Prodromal_Quality_Data[, sapply(Prodromal_Quality_Data, is.numeric)], use = "pairwise.complete.obs")
high_corr <- which(abs(cor_matrix) > 0.9 & abs(cor_matrix) < 1, arr.ind = TRUE)
data.frame(
  Var1 = rownames(cor_matrix)[high_corr[, 1]],
  Var2 = colnames(cor_matrix)[high_corr[, 2]],
  Correlation = cor_matrix[high_corr]
)
#Obtain the mean of highly correlated variables (they should be measured along the same units)
Prodromal_Quality_Data$PCT_COMBINED <- rowMeans(
  Prodromal_Quality_Data[, c("PF_READS_ALIGNED", "PF_ALIGNED_BASES")], 
  na.rm = TRUE
)

#Creating model for differential gene expression before accounting for cell value as a covariate
design_1= model.matrix(~ post_prob_v1 + gender.y + age + PCT_COMBINED + PCT_ADAPTER + MEAN_ALIGNED_READ_LENGTH, data=Prodromal_Quality_Data)
  #PCT_reads aligned and aligend bases are perfectly so I included their mean score as the the variable. 

dge <- DGEList(counts = GeneCount)
keep= filterByExpr(dge, design_1) #Filter genes based on model design
dge= dge[keep, , keep.lib.sizes=FALSE]
dge=calcNormFactors(dge) 
v=voom(dge, design_5)
fit=lmFit(v, design_5)
fit=eBayes(fit, robust=TRUE)
results=topTable(fit, coef="post_prob_v1", adjust= "BH", number =Inf) 
  #Results saved for Plotting and Gene Set Enrichment Analysis

#Identify number of genes with significant differential gene expression
sigGen= results[results$adj.P.Val<0.05, ]  

#Checking for Assumption when Testing
plot(fit$df.residual, main = "Residuals for Differential Expression") #Points are fitted on the line- assumption checked
residuals= fit$df.residual 
shapiro.test(residuals)
qqnorm(residuals)
qqline(residuals, col="red") 
#Independence between samples- Assumption checked
#Voom accounts for heteroscedasticity and it is okay if normality is deviated a bit
#BH corrections applied for multiple testing effects

#Creating model for differential gene expression after accounting for cell value as a covariate
CellFrac_DF= Cibersortx_Output_AbolsuteCellValue #Upload AbolsuteCellValue Output form Cibersortx
#Remove cell types with zero variance and/or only zeros absolute cell value across all samples
CellFrac_DF= CellFrac_DF[, -c(3, 9, 14)] 
CellFrac_DF= CellFrac_DF[, -c(14, 15, 16, 17)]

#Combining CellFrac_DF data with the Prodromal_Quality_Data data
Prodromal_Quality_Data= rownames_to_column(Prodromal_Quality_Data, var="Mixture") #Mixture refers to sample name
Prodromal_Quality_Data_all= merge(Prodromal_Quality_Data, CellFrac_DF, by="Mixture", all.x=TRUE)
rownames(Prodromal_Quality_Data_all)= Prodromal_Quality_Data_all[, 1]

GeneCount <- GeneCount[ , match(rownames(Prodromal_Quality_Data_all), colnames(GeneCount))] #Ensure the sample names remain matched between the GeneCount data and Prodromal_Quality_Data_all data
all(colnames(GeneCount) == rownames(Prodromal_Quality_Data_all))
    #Saved the Prodromal_Quality_Data_all as file

#Plotting distirbution of immune cell types 
CellFrac_DF_long= CellFrac_DF %>%
  gather(key="CellType", value="value", -Mixture)
CellFrac_DF_long$value= abs(CellFrac_DF_long$value)
ggplot(CellFrac_DF_long, aes(x=CellType, y=value))+
  geom_boxplot(aes(fill=CellType), outlier.size=2, outlier.colour="red")+
  theme_minimal()+
  labs(title="Immune Cell Types: Absolute Values ",
       x= "Cell Type",
       y= "Absolute Value")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  scale_y_continuous(limits= c(0,3.75), breaks= seq(0, 3.75, by=0.25))

#Perform correlation analyses between each cell type value and prodromal PD probability 
p_values= c()
correlations= c()
cell_type= c("T.effector.CD4.CD8_T.memory.CD8", "T.naive.CD4.CD8_T.memory.CD4", "Neutrophils.mature", "Monocytes.progenitor", "Monocytes.non.classical", "Monocytes.classical", "Erythrocytes", "NK.cells", "NK_B.progenitor", "B.cells.regulatory", "Plasma.cells", "B.cells.naive")
cell_type_2= c("Lymp_sig","T.effector.CD4.CD8_T.memory.CD8", "T.naive.CD4.CD8_T.memory.CD4", "Neutrophils.mature", "Monocytes.progenitor", "Monocytes.non.classical", "Monocytes.classical", "Erythrocytes", "NK.cells", "NK_B.progenitor", "B.cells.regulatory", "Plasma.cells", "B.cells.naive")
  #Run the below using cell_type_1 and then using cell_type_2
for (cell in cell_type_1) {
  cor_test= cor.test(abs(Prodromal_Quality_Data_all[[cell]]), Prodromal_Quality_Data_all$post_prob_v1, method="kendall")
  correlations= c(correlations, cor_test$estimate)
  p_values= c(p_values, cor_test$p.value)
}

adjusted_p= p.adjust(p_values, method="BH")
correlations_res= data.frame(
  CellType= cell_type,
  Correlation= correlations,
  PValue= p_values,
  AdjustedPValue= adjusted_p
)

print(correlations_res)

ggplot(correlations_res, aes(x=reorder(CellType, Correlation), y= Correlation, fill= AdjustedPValue<0.05))+
  geom_bar(stat= "identity", show.legend=FALSE)+
  coord_flip()+
  labs(title="Kendall Tau Correlations Between Cell Type Absolute Values and Prodromal Risk Score",
       x= "Cell Type",
       y= "Correlation Value")+
  scale_fill_manual(values= c("grey", "blue"))+
  theme_minimal()+
  theme(
    axis.text.x= element_text(angle=45, hjust=1),
    legend.title= element_text(size=17, face="bold"),
    legend.text= element_text(size=15, face="bold"))


Prodromal_Quality_Data_all= Prodromal_Quality_Data_all[, -1]

Prodromal_Quality_Data_all$post_prob_v1= as.numeric(Prodromal_Quality_Data_all$post_prob_v1)
Prodromal_Quality_Data_all$post_prob_v2= as.numeric(Prodromal_Quality_Data_all$post_prob_v2)
Prodromal_Quality_Data_all$age= as.numeric(Prodromal_Quality_Data_all$age)
Prodromal_Quality_Data_all$gender.y= factor(Prodromal_Quality_Data_all$gender.y, levels = c("male", "female"))
Prodromal_Quality_Data_all$T.naive.CD4.CD8_T.memory.CD4= as.numeric(Prodromal_Quality_Data_all$T.naive.CD4.CD8_T.memory.CD4)

Prodromal_Quality_Data_all$PF_ALIGNED_BASES= as.numeric(Prodromal_Quality_Data_all$PF_ALIGNED_BASES)
Prodromal_Quality_Data_all$PCT_PF_READS_ALIGNED= as.numeric(Prodromal_Quality_Data_all$PF_ALIGNED_BASES)
Prodromal_Quality_Data_all$PCT_ADAPTER= as.numeric(Prodromal_Quality_Data_all$PCT_ADAPTER)
Prodromal_Quality_Data_all$MEAN_ALIGNED_READ_LENGTH= as.numeric(Prodromal_Quality_Data_all$MEAN_ALIGNED_READ_LENGTH)


vars=c("gender.y", "age" , "PF_READS_ALIGNED" , "PF_ALIGNED_BASES" , "PCT_ADAPTER" , "MEAN_ALIGNED_READ_LENGTH" , "post_prob_v1", "T.naive.CD4.CD8_T.memory.CD4")

Prodromal_Quality_Data= Prodromal_Quality_Data_all[ , colnames(Prodromal_Quality_Data_all) %in% vars]

Prodromal_Quality_Data$PCT_COMBINED <- rowMeans(
  Prodromal_Quality_Data_5[, c("PF_READS_ALIGNED", "PF_ALIGNED_BASES")], 
  na.rm = TRUE
)

design_2= model.matrix(~ post_prob_v1 + gender.y + age + PCT_COMBINED + PCT_ADAPTER + MEAN_ALIGNED_READ_LENGTH + T.naive.CD4.CD8_T.memory.CD4, data=Prodromal_Quality_Data)

dge <- DGEList(counts = GeneCount)
keep= filterByExpr(dge, design_2)
dge= dge[keep, , keep.lib.sizes=FALSE]
dge=calcNormFactors(dge)
v=voom(dge, design_5)
fit=lmFit(v, design_5)
fit=eBayes(fit, robust=TRUE)
results=topTable(fit, coef="post_prob_v1", adjust= "BH", number =Inf) #Saved this file for gene set enrichment analyses
sigGen= results[results$adj.P.Val<0.05, ] 

