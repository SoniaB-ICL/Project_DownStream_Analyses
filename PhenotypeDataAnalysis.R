#Load Libraries
install.packages("patchwork")
library(patchwork)
library(ggplot2)
install.packages(c("ggplot2", "scales", "ggthemes"))
install.packages("ggstats")
library(scales)
library(ggthemes)
library(tidyr)
library(dplyr)
library(ggpubr)
install.packages("ggExtra")
library(ggExtra)
library(corrplot)

meta = read.csv("pathway/Phenotype_QC_UP.csv")
meta= meta[, -1]

meta= meta %>%
  filter(study_arm !="PD")

#Risk scores are NOT normal- so must use non-parametric testing
shapiro.test(meta$post_prob_v1)
shapiro.test(meta$post_prob_v2)

#Used Kendall instead of spearman because there are ties in risk score between samples 
corr_Risks_test= cor.test(meta$post_prob_v1, meta$post_prob_v2, method= "kendall")
print(corr_Risks_test)

#Visualise variance and distribution of risk score data
ggplot(meta, aes(x=" ", y= post_prob_v1 )) +
  geom_boxplot(color="blue", fill="lightblue") +
  labs(title= "Boxplot of Risk Score", x="", y="post_prob_v1") +
  theme_minimal()

ggplot(meta, aes(x=" ", y= post_prob_v2 )) +
  geom_boxplot(color="blue", fill="lightblue") +
  labs(title= "Boxplot of Risk Score", x="", y="post_prob_v2") +
  theme_minimal()
#do density plot or histogram insteda of box plot 
var(meta$post_prob_v1)

var(meta$post_prob_v2)

summary(meta$post_prob_v1)

summary(meta$post_prob_v2)



#Coefficient of Variation expresses how variables (%) data are relative to the mean
cv1= (sd(meta$post_prob_v1)/mean(meta$post_prob_v1)) * 100
print(cv1) 
cv2= (sd(meta$post_prob_v2)/mean(meta$post_prob_v2)) * 100
print(cv2) 

#Pairwise Variability
ggplot(meta, aes(x=post_prob_v1, y=post_prob_v2))+
  geom_point()+
  labs(title="Scatter plot of V1 vs V2", x= "post_prob_v1", y= "post_prob_v2") +
  theme_minimal()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Looking at heterogeneity in data 

#Age
corr_Risks= cor.test(meta$post_prob_v1, meta$age, method= "kendall")
print(corr_Risks) 

corr_Risks_test= cor.test(meta$age, meta$post_prob_v2, method= "kendall")
print(corr_Risks_test) 

var(meta$age)

sqrt(66.24425) 
summary(meta$age)

corr_Risks= cor.test(meta$post_prob_v1, meta$disease_duration_diag, method= "kendall")
print(corr_Risks) 

corr_Risks= cor.test(meta$post_prob_v2, meta$disease_duration_diag, method= "kendall")
print(corr_Risks)

var(meta$disease_duration_diag)

summary(meta$disease_duration_diag)


selected_vars= meta[, c("post_prob_v1", "post_prob_v2", "age")]

library(tidyr)

data_long= gather(selected_vars, key="Variable", value= "Value", post_prob_v1, post_prob_v2, age)
ggplot(data_long, aes(x=Value, fill=Variable))+
  geom_density(alpha=0.5)+
  facet_wrap(~Variable, scales="free")+
  ggtitle ("Density Plots for Prodromal Risk Scores and Age") +
  theme_minimal()+
  theme(legend.position="none")

#Possible subthreshold parkinsonism
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_motor, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_motor, method= "kendall")
print(corr_Risks) 

var(meta$LR_motor)
summary(meta$LR_motor)
 

#Olfactory Dysfunc.
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_smell, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_smell, method= "kendall")
print(corr_Risks) 

var(meta$LR_smell)

summary(meta$LR_smell)


#Blood pressure: symptomatic hypotension
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_hypo, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_hypo, method= "kendall")
print(corr_Risks) 

var(meta$LR_hypo)

summary(meta$LR_hypo)


#Urinary Dynsfunc.
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_urine, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_urine, method= "kendall")
print(corr_Risks)

var(meta$LR_urine)
summary(meta$LR_urine)


#Depression
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_depression, method= "kendall")
print(corr_Risks)
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_depression, method= "kendall")
print(corr_Risks)

var(meta$LR_depression)

summary(meta$LR_depression)


#Constipation
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_constip, method= "kendall")
print(corr_Risks)
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_constip, method= "kendall")
print(corr_Risks)

var(meta$LR_constip)

summary(meta$LR_constip)


#Daytime somnolence
corr_Risks= cor.test(meta$post_prob_v1, meta$LR_ess, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$LR_ess, method= "kendall")
print(corr_Risks)

var(meta$LR_ess)

summary(meta$LR_ess)


corMat= cor(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")], method="kendall")
p_matrix= matrix(NA, nrow=ncol(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")]), ncol=ncol(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")]))
rownames(p_matrix)= colnames(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")])
colnames(p_matrix)= colnames(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")])

selected_vars= metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")]
for (i in 1:ncol(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")])) {
  for (j in 1:ncol(metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")])) {
    test= cor.test(selected_vars[[i]], selected_vars[[j]], method="kendall")
    p_matrix[i, j] = test$p.value
  }
}
significance= ifelse(p_matrix<0.05, "*", "")

layout(matrix(c(1,2), nrow=2), heights = c(1,7))

par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "Clinical Marker Likelihood Ratios and Prodromal Risk Scores Correlation Matrix", cex=1.5, font=2)

par(mar=c(5, 5, 3, 1))
corrplot(corMat, method="color",
         type= "full",
         tl.col="black",
         tl.srt=45,
         tl.pos= "lt",
         cl.pos="r",
         p.mat=p_matrix,
         sig.level=0.05,
         insig="label_sig",
         col= colorRampPalette(c("red", "white", "blue"))(200),
         cex.lab=0.6,
         cex.axis= 0.6,
         tl.cex=0.8,
         #main="Clinical Marker Likelihood Ratios and Prodroma Risk Scores Correlation Matrix",
         addgrid.col= "grey")


selected_vars= metaRBD[, c("post_prob_v1", "post_prob_v2", "age", "LR_ess", "LR_constip", "LR_depression", "LR_urine", "LR_hypo", "LR_smell", "LR_motor")]

data_long= gather(selected_vars, key="Variable", value= "Value", post_prob_v1, post_prob_v2, age, LR_ess, LR_constip, LR_depression, LR_urine, LR_hypo, LR_smell, LR_motor)
ggplot(data_long, aes(x=Value, fill=Variable))+
  geom_density(alpha=0.5)+
  facet_wrap(~Variable, scales="free")+
  ggtitle ("Density Plots for Prodromal Risk Scores and Clinical Markers") +
  theme_minimal()+
  theme(legend.position="none")

#~~~~~~~~~~~~~~~~~~~~~~Factors not part of LR variables

corr_Risks= cor.test(meta$post_prob_v1, meta$UPDRS_III, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$UPDRS_III, method= "kendall")
print(corr_Risks) 

var(meta$UPDRS_III, na.rm=TRUE)

summary(meta$UPDRS_III)
 

corr_Risks= cor.test(meta$post_prob_v1, meta$purdue_total, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$purdue_total, method= "kendall")
print(corr_Risks) 

var(meta$purdue_total, na.rm=TRUE)

summary(meta$purdue_total)


corr_Risks= cor.test(meta$post_prob_v1, meta$sniffin_total, method= "kendall")
print(corr_Risks) 
corr_Risks= cor.test(meta$post_prob_v2, meta$sniffin_total, method= "kendall")
print(corr_Risks) 

var(meta$sniffin_total, na.rm=TRUE)

summary(meta$sniffin_total)



#Scatter Plot & Hitogram Plots
selected_vars= meta[, c("post_prob_v1", "age", "sniffin_total", "purdue_total", "UPDRS_III")]

selected_vars= na.omit(selected_vars)


layout(matrix(1:25, nrow = 5, ncol = 5, byrow = TRUE), widths = c(4,4,4,4,4), heights = rep(4,5))
# Axis label variables
vars <- colnames(selected_vars)

par(cex.main = 1.5,  
    cex.axis = 1.2, 
    cex.lab = 1.3)

par(mar = c(3, 3, 1, 1))  
par(oma = c(0.5, 0.5, 5, 0.5))

# Loop to go through each variable pair to create scatter plots as one figure
for (i in 1:5) {
  for (j in 1:5) {
    #Scatter plots
    if (i > j) {
      test= cor.test(selected_vars[, j], selected_vars[, i], method="kendall")
      tau= round(test$estimate, 2)
      p_val= sprintf("%.3E", test$p.value)
      plot(selected_vars[, j], selected_vars[, i], xlab = "", ylab = "", main = paste0("t=", tau, ", p=", p_val), 
           pch = 16, col = rgb(0.1, 0.2, 0.5, 0.5), cex = 1, cex.main=1)
      
    
      if (i == 5) { 
        mtext(vars[j], side = 1, line = 2, cex=0.8)
      }
      if (j == 1) {  
        mtext(vars[i], side = 2, line = 2, cex=0.8)
      }
    } else if (i == j) { 
      # Histogram for every variable, place on diagonal
      data_range <- range(selected_vars[, i], na.rm = TRUE)
      plot(density(selected_vars[, i], na.rm = TRUE), main = "",  
           xlab = vars[i], col = rgb(0.1, 0.2, 0.5, 0.5), 
           lwd = 2, 
           xlim = c(data_range[1] - 0.1 * diff(data_range),  
                    data_range[2] + 0.1 * diff(data_range)))
      if (i == 5) {  
        mtext(vars[i], side = 1, line = 2, cex = 0.8)
      }
      if (j == 1) {  #Histogram label
        mtext(vars[i], side = 2, line = 2, cex = 0.8)
      }
    } else {
      plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
    }
  }
}

mtext("Longitudinal Factor Distributions and Correlations", side = 3, line = 0, outer = TRUE, cex = 1)

# Reset the plot layout
par(mfrow = c(1, 1))



corMat= cor(selected_vars, use= "pairwise.complete.obs", method="kendall")
p_matrix= matrix(NA, nrow=ncol(selected_vars), ncol=ncol(selected_vars))
rownames(p_matrix)= colnames(selected_vars)
colnames(p_matrix)= colnames(selected_vars)

for (i in 1:ncol(selected_vars)) {
  for (j in 1:ncol(selected_vars)) {
    test= cor.test(selected_vars[[i]], selected_vars[[j]], method="kendall")
    p_matrix[i, j] = test$p.value
    #p_matrix[j, i] = test$p.value
  }
}

significance= ifelse(p_matrix<0.05, "*", "")

layout(matrix(c(1,2), nrow=2), heights = c(1,7))

par(mar=c(0, 0, 0, 0))
plot.new()
text(0.5, 0.5, "Logitudinal Factors and Prodromal Risk Score Correlation Matrix", cex=1.5, font=2)


par(mar=c(5, 5, 3, 1))
corrplot(corMat, method="color",
         type= "full",
         tl.col="black",
         tl.srt=45,
         tl.pos= "lt",
         cl.pos="r",
         p.mat=p_matrix,
         sig.level=0.05,
         insig="label_sig",
         col= colorRampPalette(c("red", "white", "blue"))(200),
         cex.lab=0.6,
         cex.axis= 0.6,
         tl.cex=0.8,
         addgrid.col= "grey")
mtext("Logitudinal Factors and Prodromal Risk Score Correlation Matrix", side=3, line=6, cex=1.5, font=2)

title(main="Logitudinal Factors and Prodromal Risk Score Correlation Matrix", side=3, line=6, cex=1.5, font=2)
#Put p-values in the heatmap and the legend is the corr coefficient 
 

library(tidyr)
data_long= gather(selected_vars, key="Variable", value= "Value", post_prob_v1, post_prob_v2, age, sniffin_total, purdue_total, UPDRS_III)
ggplot(data_long, aes(x=Value, fill=Variable))+
  geom_density(alpha=0.5)+
  facet_wrap(~Variable, scales="free")+
  ggtitle ("Density Plots for Prodromal Risk Scores and Longitudinal Factors") +
  theme_minimal()+
  theme(legend.position="none")


mean(Phenotype_QC_all$age) 
sd(Phenotype_QC_all$age) 

mean(Phenotype_QC_all$PCT_ADAPTER) 
sd(Phenotype_QC_all$PCT_ADAPTER) 