##1k_Microarray_IBD_Granlund.R
### Granlund_Crohns: http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0056818

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(lumi); library(limma); library(biomaRt); library(plyr); library(sva)

#1) Load Data
#------------
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_IBD_Granlund_normalized.RData")) {
  ## Load MetaData
  ## ---------------
  datMeta = read.delim("data/raw_data/Microarray/Granlund/E-MTAB-184.sdrf.txt")
  datMeta$Group = gsub("unaffected", "CTL", datMeta$Characteristics.DiseaseState.)
  datMeta$Group = gsub("diseased", "IBD", datMeta$Group)
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "IBD"))
  datMeta$Sex = gsub("female", "F", datMeta$Characteristics.Sex.)
  datMeta$Sex = gsub("male", "M", datMeta$Sex)
  datMeta$Sex = factor(datMeta$Sex, levels=c("M", "F"))
  rownames(datMeta) = datMeta$Sample.Name
  datMeta$Study="Granlund.IBD"

  ## Load Expression Data
  ## ---------------------
  data.lumi = lumiR("data/raw_data/Microarray/Granlund/atle_mage_ml_atle_mage_ml_nonorm_nobkgd_ArrayExpress_DataFile.txt")
  datExpr = lumiN(data.lumi, method="quantile")
  datExpr = log2(exprs(datExpr))
  idx = match(colnames(datExpr), datMeta$Sample.Name)
  datMeta = datMeta[idx,]

  ## Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",  host="feb2014.archive.ensembl.org")
  #f = listFilters(ensembl);   a = listAttributes(ensembl)
  identifier <- "illumina_humanht_12_v3"
  getinfo <- c("illumina_humanht_12_v3", "ensembl_gene_id","external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr),geneDat[,1])
  datProbes = geneDat[idx,]

  save(file="data/working_data/Microarray/01_Normalized/Microarray_IBD_Granlund_normalized.RData",datExpr,datProbes,datMeta)
}

load("data/working_data/Microarray/01_Normalized/Microarray_IBD_Granlund_normalized.RData")

# Remove singular batches, remove confounders --> none

# Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr_wOutliers = datExpr
datMeta_wOutliers = datMeta
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

## Batch Correction --> no batch information

to_keep = !is.na(datMeta$Sex)
datMeta = datMeta[to_keep,]
datExpr =datExpr[,to_keep]

save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_IBD_Granlund_normalized_balanced.RData",datExpr,datProbes,datMeta)

# Collapse Rows
realGenes = !is.na(datProbes[,2])
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes[,2], rowID = datProbes[,1])
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes[,1])
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes[,2]
dim(datExpr)

#Regress covariate
X = model.matrix(~Group + Sex-1, datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = X[,3] %*% t(beta[3,])
datExpr = datExpr - t(to_regress)

save(file="data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_IBD_Granlund_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)
