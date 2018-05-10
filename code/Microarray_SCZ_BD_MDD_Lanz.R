rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

set.seed(100)

if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_MDD_Lanz_normalized.RData")) {

  pfc_samples = list.files("data/raw_data/Microarray/Lanz_GSE53987/Lanz_GSE53987_RAW/Lanz_GSE53987_RAW/")
  pfc_samples = pfc_samples[grep("B46", pfc_samples)]
  data.affy = ReadAffy(celfile.path="data/raw_data/Microarray/Lanz_GSE53987/Lanz_GSE53987_RAW/Lanz_GSE53987_RAW/", filenames=pfc_samples)
  datExpr=exprs(data.affy)

  datMeta = read.csv("data/raw_data/Microarray/Lanz_GSE53987/Lanz_GSE53987_datMeta.csv")
  rownames(datMeta) = datMeta$Chip
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Race = as.factor(datMeta$Race)
  datMeta$Region[datMeta$BrainRegion=="hippocampus"] = "HC";   datMeta$Region[datMeta$BrainRegion=="Pre-frontal cortex (BA46)"] = "PFC";   datMeta$Region[datMeta$BrainRegion=="Associative striatum"] = "STR"
  datMeta$Region = as.factor(datMeta$Region)
  datMeta$Group[datMeta$Disorder=="bipolar disorder"]="BD"
  datMeta$Group[datMeta$Disorder=="control"]="CTL"
  datMeta$Group[datMeta$Disorder=="major depressive disorder"]="MDD"
  datMeta$Group[datMeta$Disorder=="schizophrenia"]="SCZ"
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "BD", "SCZ", "MDD"))
  datMeta$Study = "SCZ.BD.MDD.Lanz"
  datMeta$Subject = paste(datMeta$BrainBank, datMeta$Group, datMeta$Sex, datMeta$Age, datMeta$PMI, datMeta$Race, sep=".")

  samples = substr(colnames(datExpr),1,10)
  datMeta = datMeta[samples,]
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = sd

  RNAdeg = AffyRNAdeg(data.affy)
  datMeta$RNAdeg = RNAdeg$slope

  # Normalize
  datExpr = rma(data.affy, background =T, normalize=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr),1,10)

  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)

  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  save(file="data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_MDD_Lanz_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_MDD_Lanz_normalized.RData")

#Remove singular batches
table(datMeta$Batch)

#Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*cor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

#Balance groups
to_remove = (datMeta$Group=="SCZ") & (datMeta$RNAdeg > 5)
datExpr = datExpr[,!to_remove]
datMeta = datMeta[!to_remove,]

save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection//Microarray_SCZ_BD_MDD_Lanz_normalized_balanced.RData", datExpr, datMeta, datProbes)


#Batch Correction
mod = model.matrix(~datMeta$Group)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch, mod)
datExpr = datExpr.combat

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133_plus_2)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133_plus_2)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id

X = model.matrix(~Group+Sex+Age+PMI+pH+Race+RIN+RNAdeg-1, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,5:11]) %*% (as.matrix(beta[5:11,])))
datExpr = datExpr - t(to_regress)

idx=  which(datMeta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
datMeta$Group.BD = datMeta$Group.SCZ = datMeta$Group
datMeta$Group.SCZ[c(which(datMeta$Group=="BD"), idx[!subset.scz])] = NA
datMeta$Group.BD[c(which(datMeta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned//Microarray_SCZ_BD_MDD_Lanz_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)






