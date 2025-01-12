#Install and load packages
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("DOSE")
BiocManager::install("EnhancedVolcano")
BiocManager::install("speckle")
BiocManager::install("limma")
BiocManager::install("edgeR")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.ed.db")
BiocManager::install("enrichplot")
BiocManager::install("fgsea")
BiocManager::install("org.Hs.eg.db")
install.packages("Seurat")
install.packages("reshape2")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("RColorBrewer")
install.packages("clustree")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("stringr")

library(DOSE)
library(EnhancedVolcano)
library(speckle)
library(limma)
library(edgeR)
library(clusterProfiler)
library(Seurat)
library(reshape2)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(clustree)
library(enrichplot)
library(fgsea)
library(AnnotationDbi)
library(ggplot2)
library(ggrepel)
library(stringr)
library(org.Hs.eg.db)

#Create Volcano Plot to Visualise Differential Gene Expression

Results_genes= data.frame(read.table("pathway/DataSet", header= TRUE)) #DataSet = gene expression table output from performing differential gene expression

logFC_threshold <- 0  
p_value_threshold <- 0.05 

#Categorise genes that are significantly upregulated, significantly downregulated and non-significant 
Results_genes$Significance <- "Non-significant"
Results_genes$Significance[Results_genes$logFC > logFC_threshold & Results_genes$adj.P.Val < p_value_threshold] <- "Upregulated"
Results_genes$Significance[Results_genes$logFC < -logFC_threshold & Results_genes$adj.P.Val < p_value_threshold] <- "Downregulated"

color_scheme <- c("Non-significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")

#Subset for the top 10 significiantly expressed genes for their names to appear on the plot
top_genes= rownames(Results_genes[order(Results_genes$adj.P.Val), ][1:10, ])

Results_genes$Gene= rownames(Results_genes)

# Create the volcano plot
volcano_plot <- ggplot(Results_genes, aes(x = logFC, y = -log10(adj.P.Val), color = Significance)) +
  geom_point(alpha = 0.8, size = 2) +  
  scale_color_manual(values = color_scheme) +  #Sets the color scheme
  theme_minimal() +  
  labs(title = "Volcano Plot of Differentially Expressed Genes in RBD Cohort",
       x = "Log2 Fold Change",
       y = "Adjusted P-value (-log10)",
       color = "Significance") +  # Labels for axes and legend
  ylim(0, 2) +
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "dark green") +  
  geom_hline(yintercept = -log10(p_value_threshold), linetype = "dashed", color = "dark green") +  
  theme(
    axis.text.y= element_text(size=16, face="bold"),
    axis.title= element_text(size=16, face="bold"),
    plot.title = element_text(face="bold"))


print(volcano_plot) 

top_genes= Results_genes[order(Results_genes$adj.P.Val), ][1:10, ]
volcano_plot +
  geom_text_repel(data=top_genes, aes(label=Gene), size=2, color="black", fontface="bold", max.overlaps=15, box.padding=1.3, point.padding=1)


#Specifically subset the gene expression dataset for significantly upregulated genes: 
significant_genes_up <- Results_genes[Results_genes$adj.P.Val < 0.05 & Results_genes$logFC>0, ]
significant_genes_up$geneName= rownames(significant_genes_up)
significant_genes_up_name= significant_genes_up$geneName

#Specifically subset the gene expression dataset for significantly downregulated genes: 
significant_genes_dwn <- Results_genes[Results_genes$adj.P.Val < 0.05 & Results_genes$logFC<0, ]
significant_genes_dwn$geneName= rownames(significant_genes_dwn)
significant_genes_dwn_name= significant_genes_dwn$geneName

#Setting Up Gene Rank for Gene Set Enrichment Analysis:
Results_genes$rank= (Results_genes$logFC) * (-log10(Results_genes$adj.P.Val)) #This was used for the gene expressions before accounting for cell value as a covariate in the differential gene expression model- see below eventually this data output will stored in "DataCombine"
Results_genes$rank= sign(Results_genes$logFC) * (-log10(Results_genes$P.Value)) #This was used for the gene expressions both before/after accounting for cell value as a covariate in the differential gene expression model- see below eventaullay these data outputs will be stored in DataCombine2 and DataCombine3

#Organsing gene ranks into a gene list in dcreasing order of gene rank
GeneList= Results_genes$rank
names(GeneList)= rownames(Results_genes)
GeneList= sort(GeneList, decreasing = TRUE)

#Check whether any ranks are tied between genes, if so ensure 
dup_rank= duplicated(Results_genes$rank)
sum(dup_rank)

Enrichment_Output= gseGO(
  geneList= gene_list,
  ont= "CC",  #This is ran for BP (Biological Processes) and MF (Molecular Functions) as well 
  keyType= "SYMBOL",  #Gene names are considered as the gene symbols
  OrgDb= org.Hs.eg.db,   #Database of Homo Sapiens Annotations used to map the gene names to the gene ontology paths 
  pAdjustMethod = "BH", #Benjamini-Hochberg to correct p-values due to multiple testing effects 
  minGSSize = 5,        #Minimum gene set size for a gene pathway 
  maxGSSize = 500,      #Maximum gene set size for a gene pathway 
  nPermSimple= 10000,  #Number of Permutations
  eps=0, #Prevents incomplete encirhment computation
)

Enrichment_Output[Enrichment_Output@result$qvalue <0.05, ] #Excludes gene pathways if the false discovery rate is not less than 0.05

CC_Result= Enrichment_Output@result #Repeat for BP and MF as the ont in the gseGO()

#Combine CC, BP, and MF gene set enrichment analyses and subset for the top 20 suppressed pathways
GO1= read.csv("CC_Result", header=TRUE, sep=",") 
GO2= read.csv("BP_Result", header=TRUE, sep=",")
GO3= read.csv("MF_Result", header=TRUE, sep=",")
DataCombine= rbind(GO1, GO2, GO3)

DataCombine=DataCombine[order(-DataCombine$NES), ] 

DataCombine_Sub= DataCombine[228:248, ]


#Combine Gene Pathways Description and associate code
DataCombine_Sub$full_name= paste(DataCombine_Sub$Description, "(", Data$ID, ")", sep="")
DataCombine_Sub$full_name= str_wrap(DataCombine_Sub$full_name, width=80)

ggplot(DataCombine_Sub, aes(x=NES, y=reorder(full_name, NES), size=setSize, color=qvalue))+
  geom_point(alpha=1)+
  theme_minimal()+
  labs(
    title="Gene Pathway Suppression in RBD Cohort (Adj.P-Values; before incl. Cell Value covariate)",
    x= "Normalised Enrichment Score",
    y= "Gene Ontology Pathways",
    size= "Number of Genes",
    color= "q-value"
  )+
  scale_color_gradient(
    low="blue",
    high="red",
    name="q-value"
  )+
  #facet_grid(.~sign)+
  scale_x_continuous(limits= c(-3, 0))+
  scale_size_continuous(range= c(4, 10))+
  #scale_color_manual(values= c("Activated"="#40E0D0", "Suppressed"= "salmon"))+
  theme(
    axis.text.y= element_text(size=16, face="bold", margin=margin(t=0, r=10, b=0, l=0)),
    axis.title= element_text(size=16, face="bold"),
    strip.text = element_text(size=16, face="bold"),
    panel.spacing.x =unit(3, "lines") )


#Categorising DataCombine2 and DataCombine3 and the combining these datatset into DataCombine4
DataCombine2= DataCombine2%>%
  mutate(dataset="Before Accounting for Cell Value")

DataCombine3= DataCombine3%>%
  mutate(dataset="After Accounting for Cell Value")

DataCombine4= bind_rows(DataCombine2, DataCombine3)
DataCombine4$full_name= paste(DataCombine4$Description, "(", DataCombine4$ID, ")", sep="")

#DataCombine4 is subset to include pathways overlapping in DataCombine
pathA= DataCombine_Sub$full_name
pathB= DataCombine4$full_name

pathwaysOverlap= intersect(pathA, pathB)

DataCombine4_sub= DataCombine4[DataCombine4$full_name %in% pathwaysOverlap, ]

ggplot(DataCombine4_sub, aes(x=NES, y=Description, color=qvalue, size=setSize))+
  geom_point(aes(shape=dataset), alpha=0.7)+ #Provide different shapes for the before gene expression values vs. the after gene expression values for including cell value as a covariate
  scale_shape_manual(values=c("Before Accounting for Cell Value"=1, "After Accounting for Cell Value"=16))+
  scale_color_gradient(low="blue", high="red")+
  labs(
    title="Gene Pathway Suppression in RBD Cohort (nominal P-values; before/after incl. Cell Value covariate)",
    x= "Normalised Enrichment Score",
    y= "Pathway Description",
    color= "q-value",
    size="Gene Number",
    shape="Dataset"
  )+
  theme_minimal()+
  theme(
    axis.text.y= element_text(size=8, face="bold"),
    legend.position = "right",
    axis.title= element_text(size=16, face="bold")
  )

