##1i) Microarray_SCZ_narayan
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

#1) Load Data
#------------
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_SCZ_Narayan_normalized.RData")) {

  #Download and ProcesPhenotypic Data
  datMeta = read.csv("data/raw_data/Microarray/Narayan_GSE21138/Narayan_GSE21138_datMeta.csv")
  rownames(datMeta) = datMeta$GE

  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "SCZ"))
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Duration = as.numeric(datMeta$Duration)
  datMeta$Study = "Narayan_SCZ"
  datMeta$Subject = paste("VictorianBrainBank", datMeta$Group, datMeta$Age, datMeta$Sex,datMeta$PMI, sep="_")

  #Read in Raw Expression Data
  data.affy = ReadAffy(celfile.path="data/raw_data/Microarray/Narayan_GSE21138/Narayan_GSE21138_RAW/")

  RNAdeg = AffyRNAdeg(data.affy)
  datMeta$RNAdeg = RNAdeg$slope
  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)

  ## RMA Normalize
  datExpr = rma(data.affy, background =T, normalize=T, verbose=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substr(colnames(datExpr),1,9)
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]

  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133_plus_2"
  getinfo <- c("affy_hg_u133_plus_2", "ensembl_gene_id", "entrezgene", "external_gene_id","chromosome_name","start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)

  idx = match(rownames(datExpr), geneDat$affy_hg_u133_plus_2)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  save(file="data/working_data/Microarray/01_Normalized/Microarray_SCZ_Narayan_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_SCZ_Narayan_normalized.RData")

#Balance Groups by covariates, remove singular batches
table(datMeta$Batch)
to_keep = (datMeta$Batch != "06/14/06") & (datMeta$RNAdeg < 5.5)
datMeta = datMeta[to_keep, ]; datExpr = datExpr[,to_keep]
datMeta$Batch = factor(datMeta$Batch)

##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_Narayan_normalized_balanced.RData", datExpr, datMeta, datProbes)

## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
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

#Regress covariates
X = model.matrix(~Group+Sex+Age+PMI+pH+RNAdeg-1, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))  # Technical Covariates
datExpr = datExpr - t(to_regress)

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_Narayan_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

