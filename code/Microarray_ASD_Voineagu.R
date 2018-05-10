library(WGCNA); library(lumi); library(ggplot2); library(nlme)
library(gplots); library(GEOquery)
library(biomaRt); library(sva)

# Step 1) Download and normalize raw microarray data
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_ASD_voineagu_normalized.RData")){
  #-----------Load Raw Data----------------------------------------
  if(!file.exists("data/raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_non-normalized_RawData.csv")){
    download.file("https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE28521&format=file&file=GSE28521%5Fnon%2Dnormalized%5Fdata%2Etxt%2Egz",
                  destfile="data/raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_non-normalized_RawData.csv")}

  data.lumi = lumiR("data/raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_non-normalized_RawData.csv")
  datMeta = read.csv("data/raw_data/Microarray/Voineagu_GSE28521/Voineagu_GSE28521_metaData.csv")
  matchSN = match(sampleNames(data.lumi), datMeta$GEO_SampleName)
  datMeta = datMeta[matchSN,]

  # log2 transform and create separate datasets for: all samples, cortex, cerebellum
  dataAll.lumi<-lumiT(data.lumi, method="log2");
  dataCTX.lumi<-dataAll.lumi[,datMeta$Brain.area!="C"];
  dataCBL.lumi<-dataAll.lumi[,datMeta$Brain.area=="C"]

  #Normalize
  dataAll_N.lumi<-lumiN(dataAll.lumi, method="quantile");
  dataCTX_N.lumi<- lumiN(dataCTX.lumi, method="quantile");
  dataCBL_N.lumi<- lumiN(dataCBL.lumi, method="quantile");

  #Extract expression data for Cortex Only
  datExpr = exprs(dataCTX_N.lumi)
  datMeta = datMeta[datMeta$Brain.area!="C",]
  datExpr.prenorm = exprs(dataCTX.lumi)

  ##Re-annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "illumina_humanref_8_v3"
  getinfo <- c("illumina_humanref_8_v3", "ensembl_gene_id","external_gene_id", "entrezgene", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat[,"illumina_humanref_8_v3"])
  datProbes = cbind(rownames(datExpr), geneDat[idx,])
  rownames(datProbes) = datProbes[,1]

  #Clean Meta-Data
  datMeta$PMI[is.na(datMeta$PMI)]=mean(datMeta$PMI, na.rm=T)
  datMeta$A.C = factor(datMeta$A.C, levels = c("C", "A"))
  datMeta$Brain.area = factor(datMeta$Brain.area, levels = c("F", "T"))
  datMeta$Chip = factor(datMeta$Chip)
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Group[datMeta$A.C=="A"] = "ASD"
  datMeta$Group[datMeta$A.C=="C"] = "CTL"
  datMeta$Group = factor(datMeta$Group, levels=c("CTL", "ASD"))
  datMeta$Study = "ASD.voineagu"

  save(file="data/working_data/Microarray/01_Normalized/Microarray_ASD_voineagu_normalized.RData", datExpr, datMeta, datProbes)
}

# Step 2) Balance, outlier removal, batch correction
load("data/working_data/Microarray/01_Normalized/Microarray_ASD_voineagu_normalized.RData")

#Remove singular batches
table(datMeta$Chip)
to_remove = (datMeta$Chip == "4936551002")
datExpr = datExpr[,!to_remove]
datMeta = datMeta[!to_remove,]

#Re-balance groups
to_keep = datMeta$Age > 5 & datMeta$Age < 55
datMeta = datMeta[to_keep,]
datExpr = datExpr[,to_keep]

datMeta$Chip = factor(datMeta$Chip)

## Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$Batch = datMeta$Chip
save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_ASD_voineagu_normalized_balanced.RData", datExpr, datMeta, datProbes)

## Batch Correction
mod = model.matrix(~A.C, data=datMeta)
batch = factor(datMeta$Chip)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat

## Collapse Probes to Genes
realGenes = !is.na(datProbes$ensembl_gene_id)  #
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$illumina_humanref_8_v3)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$illumina_humanref_8_v3)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

## Regress all Technical and Non-DX biological covariates
rownames(datMeta) = datMeta$GEO_GSM
colnames(datExpr) = rownames(datMeta)

X = model.matrix(~Group+Brain.area+Sex+Age+RIN+PMI-1, data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = t(as.matrix(X[,3:7]) %*% (as.matrix(beta[3:7,])))
datExpr = datExpr - to_regress

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_voineagu_normalized_CR_cleaned.RData", datExpr, datMeta, datProbes)






