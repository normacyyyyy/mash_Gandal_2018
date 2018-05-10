##1l) Microarray_IBD_Noble
#--------------------

rm(list=ls())
options(stringsAsFactors=FALSE)
#source("http://bioconductor.org/biocLite.R")
library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

#1) Load Data
#------------
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData")) {

  datMeta=read.csv("data/raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_datMeta2.csv")
  rownames(datMeta)= datMeta$GSM
  datMeta$Sex = as.factor(datMeta$Gender)
  datMeta$Ethnicity = as.factor(datMeta$Ethnicity)
  datMeta$Group = as.factor(gsub("Normal", "CTL", gsub("UC", "IBD", datMeta$Disease)))
  datMeta$Tissue = as.factor(datMeta$Anatomic_Location)
  datMeta$Inflammation_State = as.factor(datMeta$Inflammation_State)
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$Study= "Noble.IBD"
  datMeta$Array = "Agilent_G4112A"
  datMeta$Batch = as.factor(datMeta$Run_Date)

  ## Read in Expression Data
  filenames=list.files("data/raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_RAW/")
  RG = read.maimages(files=filenames,source="agilent",path="data/raw_data/Microarray/Noble_GSE11223/Noble_GSE11223_RAW")
  RGb = backgroundCorrect(RG,method= "normexp",offset=50)
  MA <- normalizeWithinArrays(RGb,method="loess")
  MA.q <- normalizeBetweenArrays(MA,method="quantile")
  datExpr = getEAWP(MA.q)$exprs
  rownames(datExpr) = getEAWP(MA.q)$probes$ProbeName
  colnames(datExpr) = gsub(".txt", "", colnames(datExpr))

  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]

  to_remove = which(datMeta$Gender=="unknown")
  datMeta = datMeta[-to_remove,]
  datExpr = datExpr[,-to_remove]


  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  a = listAttributes(ensembl); f=listFilters(ensembl)
  identifier <- "efg_agilent_wholegenome_4x44k_v1"
  getinfo <- c("efg_agilent_wholegenome_4x44k_v1", "ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$efg_agilent_wholegenome_4x44k_v1)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  save(file="data/working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_IBD_noble_normalized.RData")

## Remove singular batch --> no singular batches
table(datMeta$Batch)

## Remove group x covariate confounds
to_keep = (datMeta$Ethnicity != "ASIAN") & (datMeta$Ethnicity != "JEWISH")

to_keep = to_keep & !((datMeta$Group=="IBD") & (datMeta$Sex=="M") & (datMeta$Age> 50))
to_keep = to_keep & !((datMeta$Group=="CTL") & (datMeta$Sex=="F") & (datMeta$Age< 30))
table(to_keep)
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]
datMeta$Ethnicity = factor(datMeta$Ethnicity)
datMeta$Sex = factor(datMeta$Sex)

## Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_IBD_noble_normalized_balanced.RData", datExpr, datMeta, datProbes)


## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr = ComBat(datExpr, batch=batch, mod=mod)

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  & !duplicated(rownames(datExpr))
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$efg_agilent_wholegenome_4x44k_v1)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$efg_agilent_wholegenome_4x44k_v1)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id


# Regress Covariates
X = model.matrix(~Group + Sex + Age + Tissue-1, datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))
datExpr = datExpr - t(to_regress)

save(file="data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_IBD_noble_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

