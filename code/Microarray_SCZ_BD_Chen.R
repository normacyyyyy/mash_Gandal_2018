library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

set.seed(100)
if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData")) {

  datMeta=read.csv("data/raw_data/Microarray/Chen_GSE35978/GSE35978_datMeta.csv")
  rownames(datMeta)= datMeta$Chip
  datMeta$Region = as.factor(datMeta$Region)
  datMeta$Group[datMeta$Group=="Bipolar"] = "BD"; datMeta$Group[datMeta$Group=="Control"] = "CTL";  datMeta$Group[datMeta$Group=="Schizophrenia"] = "SCZ";   datMeta$Group[datMeta$Group=="Depression"] = "MDD";
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "MDD", "BD", "SCZ"))
  datMeta$Age = as.numeric(datMeta$Age)
  datMeta$PMI = as.numeric(datMeta$PMI)
  datMeta$pH = as.numeric(datMeta$pH)
  datMeta$Sex=as.factor(datMeta$Sex)
  datMeta$Study="SCZ.BD.Chen"

  #Remove NA group
  to_keep = !is.na(datMeta$Group)
  datMeta = datMeta[to_keep,]

  ## Get parietal cortex samples only
  ctx_only = rownames(datMeta)[datMeta$Region=="PCTX"]
  all_samples = list.files("data/raw_data/Microarray/Chen_GSE35978/GSE35978_CEL/")
  ctx_sample_files = all_samples[pmatch(ctx_only, all_samples)]

  ## Read in Expression Data
  data.affy = ReadAffy(celfile.path="data/raw_data/Microarray/Chen_GSE35978/GSE35978_CEL", filenames=ctx_sample_files)
  datExpr = affy::rma(data.affy,verbose=T,normalize=T,background=T)
  datExpr = exprs(datExpr)
  colnames(datExpr) = substring(colnames(datExpr), 1, 9)

  ## Get RNA degradation
  RNAdeg = AffyRNAdeg(data.affy)
  RNAdeg$sample.names = substring(RNAdeg$sample.names,1,9)
  idx=  match(rownames(datMeta), RNAdeg$sample.names)
  datMeta$RNAdeg = RNAdeg$slope[idx]

  ## Get batch information
  batch = as.factor(substring(protocolData(data.affy)$ScanDate,1,10))
  idx = match(rownames(datMeta), colnames(datExpr))
  datMeta$Batch = batch[idx]

  ## Align expression and phenoData matricies
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]

  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  a = listAttributes(ensembl)
  identifier <- "affy_hugene_1_0_st_v1"
  getinfo <- c("affy_hugene_1_0_st_v1", "ensembl_gene_id","hgnc_symbol", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$affy_hugene_1_0_st_v1)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  save(file="data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_chen_normalized.RData")

## Remove singular batch
table(datMeta$Batch)
to_keep = datMeta$Batch != "2010-04-14"
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]

## Remove MDD group (confounded by pH)
to_keep = datMeta$Group != "MDD"
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]

## Remove samples to match PMI, pH
to_keep = datMeta$Group!="MDD" & datMeta$PMI < 70
to_keep = to_keep & ((datMeta$pH < 6.9) | !(datMeta$Group=="CTL")) & (datMeta$pH > 5.8)
table(to_keep)
datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]
datMeta$Group = factor(datMeta$Group)

## Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$Batch = factor(datMeta$Batch)
save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_BD_chen_normalized_balanced.RData", datExpr, datMeta, datProbes)

## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr = ComBat(datExpr, batch=batch, mod=mod)

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hugene_1_0_st_v1)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hugene_1_0_st_v1)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id

# Regress Covariates
X = model.matrix(~Group +Sex + Age + pH + PMI + RNAdeg - 1, datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:8]) %*% (as.matrix(beta[4:8,])))  # PMI + Sex +  Age + RNAdeg
datExpr = datExpr - t(to_regress)

idx=  which(datMeta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
datMeta$Group.BD = datMeta$Group.SCZ = datMeta$Group
datMeta$Group.SCZ[c(which(datMeta$Group=="BD"), idx[!subset.scz])] = NA
datMeta$Group.BD[c(which(datMeta$Group=="SCZ"), idx[subset.scz])] = NA

save(file="data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_BD_Chen_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)
