library(WGCNA); library(affy); library(limma); library(biomaRt); library(sva)

if(!file.exists("data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_Iwa_normalized.RData")) {

  #Load Meta Data
  datMeta = read.csv(file="data/raw_data/Microarray/Iwa_GSE12649/Iwa_GSE12649_datMeta.csv", head=T)
  rownames(datMeta) = paste(datMeta$Filename,".CEL.gz",sep="")
  datMeta$Sex = as.factor(datMeta$Sex)
  datMeta$Race = as.factor(datMeta$Race)
  datMeta$Group = datMeta$Profile
  datMeta$Group[datMeta$Group=="Bipolar"]="BD"
  datMeta$Group[datMeta$Group=="Schizophrenia"]="SCZ"
  datMeta$Group[datMeta$Group=="Control"]="CTL"
  datMeta$Group[datMeta$Group=="Depression"]="MDD"
  datMeta$Group = factor(datMeta$Group, levels = c("CTL", "BD", "SCZ"))

  datMeta$RD = as.factor(datMeta$RD)
  datMeta = datMeta[!is.na(datMeta$Group),]
  datMeta$Study="SCZ.BD.Iwa"

  #Load Expression Data
  data.affy = ReadAffy(celfile.path="data/raw_data/Microarray/Iwa_GSE12649/Iwa_GSE12649_raw/Iwa_GSE12649_raw/", filenames = paste(datMeta$Filename,".CEL.gz",sep=""))
  datExpr = affy::rma(data.affy, background =T, normalize=T)
  datExpr = exprs(datExpr)
  idx = match(colnames(datExpr), rownames(datMeta))
  datMeta = datMeta[idx,]


  RNAdeg = AffyRNAdeg(data.affy)
  datMeta$RNAdeg = RNAdeg$slope

  sd = protocolData(data.affy)$ScanDate; sd = substring(sd,1,nchar(sd)-9)
  datMeta$Batch = as.factor(sd)
  table(sd)

  #Annotate Probes
  ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
  identifier <- "affy_hg_u133a"
  getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene", "external_gene_id", "chromosome_name", "start_position", "end_position")
  geneDat <- getBM(attributes = getinfo,filters=identifier,values=rownames(datExpr),mart=ensembl)
  idx = match(rownames(datExpr), geneDat$affy_hg_u133a)
  datProbes = cbind(rownames(datExpr), geneDat[idx,])

  save(file="data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_Iwa_normalized.RData", datExpr, datMeta, datProbes)
}

load("data/working_data/Microarray/01_Normalized/Microarray_SCZ_BD_Iwa_normalized.RData")

#Balance Groups by batch, RNAdeg
table(datMeta$Batch)
to_keep = !(datMeta$Batch == "07/09/03") & (datMeta$RNAdeg < 3.8) & (datMeta$RNAdeg > 2) & !is.na(datMeta$Group) & (datMeta$Filename != "77-HGU133A")
table(to_keep)
datMeta = datMeta[to_keep,]; datExpr = datExpr[,to_keep]

##Outlier Removal
sdout <- 2; normadj <- (0.5+0.5*bicor(datExpr, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj); ku <- netsummary$Connectivity; z.ku <- (ku-mean(ku))/sqrt(var(ku))
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(datExpr)[outliers]); print(table(outliers))
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]

datMeta$Batch=factor(datMeta$Batch)
save(file="data/working_data/Microarray/02_NormalizedBalanced_noBatchCorrection/Microarray_SCZ_BD_Iwa_normalized_balanced.RData", datExpr, datMeta, datProbes)

## Batch Correction
mod = model.matrix(~Group, data=datMeta)
batch = factor(datMeta$Batch)
datExpr.combat = ComBat(datExpr, batch=batch, mod=mod)
datExpr = datExpr.combat

# Collapse Rows
realGenes = !is.na(datProbes$ensembl_gene_id)  #
table(realGenes)
datExpr = datExpr[realGenes,]
datProbes = datProbes[realGenes,]

CR = collapseRows(datExpr, rowGroup = datProbes$ensembl_gene_id, rowID = datProbes$affy_hg_u133a)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], datProbes$affy_hg_u133a)
datProbes = datProbes[idx,]
rownames(datProbes) = datProbes$ensembl_gene_id
dim(datExpr)

X = model.matrix(~Group+Sex+Age+pH+PMI+RNAdeg-1,data=datMeta)
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)

to_regress = (as.matrix(X[,4:8]) %*% (as.matrix(beta[4:8,])))
datExpr = datExpr - t(to_regress)

idx=  which(datMeta$Group=="CTL")
subset.scz = runif(length(idx)) > 0.5
datMeta$Group.BD = datMeta$Group.SCZ = datMeta$Group
datMeta$Group.SCZ[c(which(datMeta$Group=="BD"), idx[!subset.scz])] = NA
datMeta$Group.BD[c(which(datMeta$Group=="SCZ"), idx[subset.scz])] = NA

save(file = "data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_SCZ_BD_Iwa_normalized_CR_regressed.RData", datExpr, datMeta, datProbes)

