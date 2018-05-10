library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(gplots); library(GEOquery)
library(biomaRt); library(sva)

# Step 1) Download and normalize raw microarray data
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_ASD_Chow_normalized.RData")){
  #Load & clean metaData
  datMeta = read.csv("data/raw_data/Microarray/Chow_GSE28475/Chow_GSE28475_datMeta.csv")
  rownames(datMeta) = datMeta$GSM
  datMeta$Group = factor(datMeta$Group,levels =  c("CTL", "ASD"))
  datMeta$Sex = as.factor(datMeta$SEX)
  datMeta$Batch = as.factor(datMeta$Batch)
  datMeta$Age = as.numeric(datMeta$AGE)
  datMeta$ChipID = as.factor(datMeta$ChipID)
  datMeta$ChipPosition = as.factor(datMeta$ChipPosition)
  datMeta$ID = gsub(" ", ".", datMeta$Sample2)
  datMeta$Study="ASD.chow"

  datExpr = read.delim(file="data/raw_data/Microarray/Chow_GSE28475/Chow_GSE28475_quantile_normalized.txt")
  rownames(datExpr) = datExpr[,1]
  datExpr=datExpr[,-1]
  datExpr= datExpr[,seq(1,65,by=2)]
  colnames(datExpr) = gsub("DASL_Frozen_","", colnames(datExpr))

  idx = match(colnames(datExpr), datMeta$ID)
  datMeta = datMeta[idx,]

  ##Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "illumina_humanref_8_v3"
  getinfo <- c("illumina_humanref_8_v3", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)

  idx = match(rownames(datExpr), geneDat$illumina_humanref_8_v3)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  colnames(datExpr) = rownames(datMeta)
  save(file="data/working_data/Microarray/01_Normalized/Microarray_ASD_Chow_normalized.RData", datMeta, datExpr, datProbes)
}

# Step 2) Balance, outlier removal, batch correction
load("data/working_data/Microarray/01_Normalized/Microarray_ASD_Chow_normalized.RData")

sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
plot(Z.K, col = as.numeric(datMeta$Group), pch=19, main="Outlier detection", ylab="Network connectivity (z score)")
legend("bottomleft",pch=16, legend = levels(datMeta$Group), col = 1:2)
abline(h=-2, lty=2)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)

table(datMeta$Batch)
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat

save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_ASD_Chow_normalized_balanced.RData", datMeta, datExpr, datProbes)

# CollapseRows Probes --> Genes
realGenes = !is.na(datProbes$ensembl_gene_id)
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$illumina_humanref_8_v3)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$illumina_humanref_8_v3)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

## Regress all Technical and Non-DX biological covariates
datMeta$RIN[is.na(datMeta$RIN)]= mean(datMeta$RIN,na.rm=T)
X = model.matrix(~Group+Age+PMI+RIN-1, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,3:5]) %*% (as.matrix(beta[3:5,])))
datExpr = datExpr - t(to_regress)

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_chow_normalized_CR_cleaned.RData", datExpr, datMeta, datProbes)


