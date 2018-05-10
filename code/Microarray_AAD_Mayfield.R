##1h) Microarray_AAD_Mayfield
#---------------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva); library(lumi)

#1) Load Data
#------------
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_AAD_mayfield_normalized.RData")) {

  ## Format MetaData
  ## ---------------
  datMeta = read.csv("data/raw_data/Microarray/Mayfield_GSE29555/Mayfield_GSE29555_datMeta.csv")
  rownames(datMeta)=datMeta$ID
  datMeta$Group = as.factor(datMeta$Group)
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$Study = "AAD.mayfield"

  ## Format ExpressionData
  ## ---------------------
  data.lumi = lumiR("data/raw_data/Microarray/Mayfield_GSE29555/Mayfield_GSE29555_non-normalized_region1.txt")
  datExpr = lumiN(data.lumi, method="quantile")
  datExpr = log2(exprs(datExpr))

  idx = match(colnames(datExpr),rownames(datMeta));
  datMeta = datMeta[idx,]

  ## Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "illumina_humanht_12_v3"
  getinfo <- c("illumina_humanht_12_v3", "ensembl_gene_id","external_gene_id", "entrezgene","chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr),geneDat[,1])
  datProbes = geneDat[idx,]

  save(file="data/working_data/Microarray/01_Normalized/Microarray_AAD_mayfield_normalized.RData",datExpr,datProbes,datMeta)
}

load("data/working_data/Microarray/01_Normalized/Microarray_AAD_mayfield_normalized.RData")

# Balance groups
#--> No singular batches or group x covariate confounds

# Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_AAD_mayfield_normalized_balanced.RData",datExpr,datProbes,datMeta)

# Collapse Rows
realGenes = !is.na(datProbes[,2])
table(realGenes)
datExpr = datExpr[realGenes,]; datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes[,2], rowID = datProbes[,1])
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes[,1])
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes[,2]
dim(datExpr)

## Regress Covariates
X = model.matrix(~Group+Age+Sex+PMI+pH+RIN-1, data = datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))  # PMI + Sex +  Age + RNAdeg
datExpr = datExpr - t(to_regress)

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AAD_mayfield_normalized_CR_regressed.RData", datExpr,datMeta, datProbes)

