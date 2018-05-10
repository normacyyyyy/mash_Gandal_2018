library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(gplots); library(GEOquery)
library(biomaRt); library(sva); library(affy)

# Step 1) Download and normalize raw microarray data
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_ASD_Garbett_normalized.RData")) {

  datMeta = read.csv("data/raw_data/Microarray/Garbett/Garbett_datMeta.csv")
  datMeta$Study = "ASD.garbett"
  datMeta$Group[datMeta$Group=="autism"] = "ASD"
  datMeta$Group[datMeta$Group=="control"] = "CTL"
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "ASD"))
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta = datMeta[ c(8,3,1,2,9,7,6,4,11,12,10,5),]
  data.affy = ReadAffy(celfile.path = "data/raw_data/Microarray/Garbett/Garbett_CEL/")

  #Calculate Batch and 5'/3' Bias Covariates
  RNAdeg = AffyRNAdeg(data.affy)
  datMeta$RNAdeg = RNAdeg$slope
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)

  ## RMA Normalize
  datExpr = affy::rma(data.affy, background =T, normalize=T, verbose=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr), 3, as.numeric(gregexpr("_", colnames(datExpr)))-2)
  datMeta$Subject2=gsub("UMB","",datMeta$Subject2)
  rownames(datMeta) = datMeta$Subject2
  idx = match(colnames(datExpr), datMeta$Subject2)
  datMeta = datMeta[idx,]


  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)

  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  save(file= "data/working_data/Microarray/01_Normalized/Microarray_ASD_Garbett_normalized.RData", datExpr, datMeta, datProbes)
}

# Step 2) Balance, outlier removal, batch correction
load(file= "data/working_data/Microarray/01_Normalized/Microarray_ASD_Garbett_normalized.RData")

#Remove Singular Batches
table(datMeta$Batch)
to_keep = (datMeta$Batch != "06/22/06")
datMeta = datMeta[to_keep,]
datExpr = datExpr[,to_keep]
datMeta$Batch = factor(datMeta$Batch)


##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))

save(file= "data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_ASD_Garbett_normalized_balanced.RData", datExpr, datMeta, datProbes)

## Batch Correction
mod = model.matrix(~datMeta$Group)
batch = factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod,prior.plots = F)
datExpr = datExpr.combat

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  #
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133_plus_2)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133_plus_2)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

##Regress covariates
X = model.matrix(~Group+Sex+Age+PMI+RIN+RNAdeg-1, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))
datExpr = datExpr - t(to_regress)

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_Garbett_normalized_CR_regressed.RData.RData", datExpr, datMeta, datProbes)

