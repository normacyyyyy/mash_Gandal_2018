rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)


#1) Load Data
#------------
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_SCZ_Maycox_normalized.RData")) {

  data.affy = ReadAffy(celfile.path="data/raw_data/Microarray/Maycox_GSE17612/Maycox_GSE17612_RAW/")
  datExpr = rma(data.affy, normalize=T, background=T, verbose=T)
  datExpr = exprs(datExpr)

  batch = as.factor(substr(protocolData(data.affy)$ScanDate,1,8))

  datMeta = read.csv("data/raw_data/Microarray/Maycox_GSE17612/Maycox_GSE17612_datMeta.csv")
  rownames(datMeta) = datMeta$Chip
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "SCZ"))
  datMeta$Sex = as.factor(datMeta$Sex)
  RNAdeg = AffyRNAdeg(data.affy)
  datMeta$RNAdegBias = RNAdeg$slope
  datMeta$Batch = batch
  datMeta$Study = "SCZ.maycox"
  datMeta$Subject = paste(datMeta$BrainBank, datMeta$ID,sep="_")

  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id",  "chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)

  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  colnames(datExpr) = substring(colnames(datExpr),1,9)
  save(file="data/working_data/Microarray/01_Normalized/Microarray_SCZ_Maycox_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_SCZ_Maycox_normalized.RData")

## Remove singular batches
table(datMeta$Batch)

##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

#Balance Groups, Remove Singular Batches
#--> No singular batches or unbalanced covariates
save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_Maycox_normalized_balanced.RData", datExpr, datMeta, datProbes)


## Batch Correction
#plot(datMeta$Group  ~ datMeta$Batch, main="Batch x Group", xlab = "Batch", ylab="") # No group by batch confound
mod = model.matrix(~Group, data=datMeta)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
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

#Regress Covariates
X = model.matrix(~Group+Sex+Age+PMI+RNAdegBias-1, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
to_regress = (as.matrix(X[,3:6]) %*% (as.matrix(beta[3:6,])))
datExpr = datExpr - t(to_regress)

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_Maycox_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)
