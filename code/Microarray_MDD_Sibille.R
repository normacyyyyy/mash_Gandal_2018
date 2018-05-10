##1j) Microarray_MDD_Sibille
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

set.seed(100)
#1) Load Data
#------------
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_MDD_Sibille1_normalized.RData")) {

  #Download and ProcesPhenotypic Data
  ## GSE54565, GSE54567, GSE54568, GSE54571, GSE54572 -- GPL570
  datMeta = read.csv("data/raw_data/Microarray/Sibille/Sibille_datMeta.csv")
  rownames(datMeta) = datMeta$GSM
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Group = factor(gsub("Control", "CTL", datMeta$GROUP), levels = c("CTL", "MDD"))
  datMeta$Region = as.factor(datMeta$REGION)
  datMeta$StudyGSE = as.factor(datMeta$GSE)
  datMeta$Study = "MDD_Sibille"
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Race = as.factor(datMeta$Race)

  #Read in Raw Expression Data
  data.affy = ReadAffy(celfile.path="data/raw_data/Microarray/Sibille/Sibille/Sibille_GSE54565_GSE54567_GSE54568_GSE54571_GSE54572_RAW/")
  idx = match(substr(colnames(exprs(data.affy)),1,10), rownames(datMeta))
  datMeta = datMeta[idx,]

  RNAdeg = AffyRNAdeg(data.affy)
  datMeta$RNAdeg = RNAdeg$slope
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)

  ## RMA Normalize
  datExpr = rma(data.affy, background =T, normalize=T, verbose=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr),1,10)
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]

  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id","chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)

  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  save(file="data/working_data/Microarray/01_Normalized/Microarray_MDD_Sibille1_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_MDD_Sibille1_normalized.RData")

## Remove singular batch
table(datMeta$Batch)
to_keep = (datMeta$Batch != "12/03/09")
datMeta = datMeta[to_keep, ]; datExpr = datExpr[,to_keep]
datMeta$Batch = factor(datMeta$Batch)

datMeta$Region = factor(datMeta$Region)

## Remove Group x covariate confound --> none
# model = model.matrix(~Group+Age+Region+Race+Sex+PMI+pH+RIN+RNAdeg+Batch, data=datMeta)
save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_MDD_Sibille1_normalized_balanced.RData", datExpr, datMeta, datProbes)

## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat

##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133_plus_2)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133_plus_2)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

datMeta$Region = factor(datMeta$Region)
X = model.matrix(~Group+Age+Region+Race+Sex+PMI+pH+RIN+RNAdeg-1, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:10]) %*% (as.matrix(beta[3:10,])))  # Technical Covariates
datExpr = datExpr - t(to_regress)

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_MDD_Sibille1_Normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

