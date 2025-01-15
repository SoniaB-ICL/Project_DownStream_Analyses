
#Setting up mixture file and signature matrix for cell deconvolution. Ensuring each file is in the correct layout and filtered by common genes. 
GeneCount = read.csv("pathway/GeneCounts_RAW_FINAL.csv", header=TRUE, row.names=1)

GeneCount$NAME=rownames(GeneCount)
GeneCount <- GeneCount %>%
  select(NAME, everything())

rownames(GeneCount)= NULL

write.table(GeneCount, "pathway/GeneCountExport.txt", sep= "\t", row.names= FALSE, col.names=FALSE, quote=FALSE) 

GeneCountExport= read.table("pathway/GeneCount.txt", header= TRUE, sep="\t", stringsAsFactors=FALSE)


non_NUM= names(GeneCountExport)[sapply(GeneCountExport, function(col) !is.numeric(col))]
print(non_NUM)

GeneCountExport$NAME=as.character(GeneCountExport$NAME)
GeneCountExport[, -1]=lapply(GeneCountExport[, -1], function(x) as.numeric(as.character(x)))
str(GeneCountExport)

any(is.na(GeneCountExport))

write.table(GeneCountExport, "pathway/GeneCountExport.txt", sep= "\t", row.names= FALSE, col.names=TRUE, quote=FALSE) 


SigMat = read.table("pathway/cibersort_input_scWB15v4_signature_matrix.txt", header=TRUE, sep="\t")

SigMat$NAME=as.character(SigMat$NAME)
SigMat[, -1]=lapply(SigMat[, -1], function(x) as.numeric(as.character(x)))
str(SigMat)

write.table(SigMat, "pathway/SigMatExport.txt", sep= "\t", row.names= FALSE, col.names=TRUE, quote=FALSE) 

summary(SigMat)
summary(GeneCountExport)

rownames(SigMat)= SigMat[, 1]
SigMat= SigMat[, -1]
rownames(GeneCountExport)= GeneCountExport[, 1]
GeneCountExport= GeneCountExport[, -1]


common_genes= intersect(rownames(SigMat), rownames(GeneCountExport)) #Filtering by common genes 

signature_filtered= SigMat[common_genes, ]
genes_filtered= GeneCountExport[common_genes, ]

rownames(genes_filtered)= genes_filtered[, 1]
genes_filtered= genes_filtered[, -1]
cpmGene= cpm(genes_filtered)
cpmGene= as.data.frame(cpmGene)

rownames(signature_filtered)= signature_filtered[ , 1]
signature_filtered= signature_filtered[, -1]
cpmSig= cpm(signature_filtered)
cpmSig= as.data.frame(cpmSig)


cpmSig= cpmSig[match(rownames(cpmGene), rownames(cpmSig)), ]

all(rownames(cpmSig)== rownames(cpmGene))

cpmSig= rownames_to_column(cpmSig, var="NAME")
cpmGene= rownames_to_column(cpmGene, var="NAME")

write.table(cpmSig, "pathway/SigMatExportUPDATED.txt", sep= "\t", row.names= FALSE, col.names=TRUE, quote=FALSE) #Signature Matrix file
write.table(cpmGene, "pathway/GeneCountExportUPDATED.txt", sep= "\t", row.names= FALSE, col.names=TRUE, quote=FALSE) #Mixture file


