rm(list=ls()); options(stringsAsFactors=F)
#source("http://bioconductor.org/biocLite.R")

library(nlme)

# ASD
ASD_multiExpr = vector(mode="list",length = 3)
load("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_voineagu_normalized_CR_cleaned.RData")
i=1; ASD_multiExpr[[i]]$datExpr = datExpr; ASD_multiExpr[[i]]$datMeta = datMeta ;ASD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)
load("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_chow_normalized_CR_cleaned.RData")
i=2; ASD_multiExpr[[i]]$datExpr = datExpr; ASD_multiExpr[[i]]$datMeta = datMeta ;ASD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)
load("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_ASD_Garbett_normalized_CR_regressed.RData.RData")
i=3; ASD_multiExpr[[i]]$datExpr = datExpr; ASD_multiExpr[[i]]$datMeta = datMeta ;ASD_multiExpr[[i]]$datProbes = datProbes; rm(datExpr,datMeta,datProbes)

genes = intersect(intersect(rownames(ASD_multiExpr[[1]]$datExpr), rownames(ASD_multiExpr[[2]]$datExpr)), rownames(ASD_multiExpr[[3]]$datExpr))
for(i in 1:3){ASD_multiExpr[[i]]$datExpr = ASD_multiExpr[[i]]$datExpr[match(genes,rownames(ASD_multiExpr[[i]]$datExpr)),]}

ASD_datExpr = cbind(ASD_multiExpr[[1]]$datExpr,ASD_multiExpr[[2]]$datExpr,ASD_multiExpr[[3]]$datExpr)
ASD_datMeta = rbind(ASD_multiExpr[[1]]$datMeta[,c("Study", "Subject","Group")],ASD_multiExpr[[2]]$datMeta[,c("Study", "Subject","Group")],ASD_multiExpr[[3]]$datMeta[,c("Study", "Subject","Group")])

# SCZ_BD_MDD
files = dir("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("SCZ",files)]
n = length(files)
SCZ_BD_MDD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  SCZ_BD_MDD_multiExpr[[i]]$datExpr = datExpr
  SCZ_BD_MDD_multiExpr[[i]]$datMeta = datMeta
  SCZ_BD_MDD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes = intersect(rownames(SCZ_BD_MDD_multiExpr[[1]]$datExpr), rownames(SCZ_BD_MDD_multiExpr[[2]]$datExpr))
for(i in 3:n){genes = intersect(genes, rownames(SCZ_BD_MDD_multiExpr[[i]]$datExpr))}
for(i in 1:n){SCZ_BD_MDD_multiExpr[[i]]$datExpr = SCZ_BD_MDD_multiExpr[[i]]$datExpr[match(genes,rownames(SCZ_BD_MDD_multiExpr[[i]]$datExpr)),]}

SCZ_BD_MDD_datExpr = data.frame(row.names = genes)
SCZ_BD_MDD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  SCZ_BD_MDD_datExpr = cbind(SCZ_BD_MDD_datExpr, SCZ_BD_MDD_multiExpr[[i]]$datExpr)
  SCZ_BD_MDD_datMeta = rbind(SCZ_BD_MDD_datMeta, SCZ_BD_MDD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group")])
}
SCZ_BD_MDD_datMeta = SCZ_BD_MDD_datMeta[-1,]
SCZ_BD_MDD_datMeta$Group = factor(SCZ_BD_MDD_datMeta$Group)

# MDD
load("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_MDD_Sibille1_Normalized_CR_regressed.RData")
MDD_datExpr = datExpr; rm(datExpr)
MDD_datMeta = datMeta[,c("Study", "Subject","Group")]; rm(datMeta)
MDD_datProbes = datProbes; rm(datProbes)

# AAD
load("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AAD_mayfield_normalized_CR_regressed.RData")
AAD_datExpr = datExpr; rm(datExpr)
AAD_datMeta = datMeta[,c("Study", "Subject","Group")]; rm(datMeta)
AAD_datProbes = datProbes; rm(datProbes)
levels(AAD_datMeta$Group)[levels(AAD_datMeta$Group)=="ETOH"] = "AAD"

# Compile
multiExpr = vector(mode="list", length=4)
multiExpr[[1]]$datExpr = ASD_datExpr;  multiExpr[[1]]$datMeta = ASD_datMeta
multiExpr[[2]]$datExpr = SCZ_BD_MDD_datExpr; multiExpr[[2]]$datMeta = SCZ_BD_MDD_datMeta
multiExpr[[3]]$datExpr = MDD_datExpr; multiExpr[[3]]$datMeta = MDD_datMeta
multiExpr[[4]]$datExpr = AAD_datExpr; multiExpr[[4]]$datMeta = AAD_datMeta
# multiExpr[[5]]$datExpr = IBD_datExpr; multiExpr[[5]]$datMeta = IBD_datMeta
names(multiExpr) = c("ASD", "SCZ_BD_MDD", "MDD", "AAD")

genes = intersect(rownames(multiExpr[[1]]$datExpr), rownames(multiExpr[[2]]$datExpr))
for(i in 3:4){genes = intersect(genes, rownames(multiExpr[[i]]$datExpr))}

for(i in 1:4){multiExpr[[i]]$datExpr = multiExpr[[i]]$datExpr[match(genes, rownames(multiExpr[[i]]$datExpr)),]}

AllExpr = data.frame(row.names = genes)
AllMeta = data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:4) {
  AllExpr = cbind(AllExpr, multiExpr[[i]]$datExpr)
  AllMeta = rbind(AllMeta, multiExpr[[i]]$datMeta)
}
AllMeta = AllMeta[-1,]
AllMeta$Group = factor(AllMeta$Group)

# GLM
Chat = matrix(NA, nrow=nrow(AllExpr), ncol=6)
SEmtx = matrix(NA, nrow=nrow(AllExpr), ncol=6)
DFmtx = matrix(NA, nrow=nrow(AllExpr), ncol=6)
for(i in 1:nrow(AllExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(AllExpr[i,])
  tryCatch({
    tSummary = summary(lme(expr ~ Group + Study - 1, data = AllMeta, random=~1|Subject, control = list(opt = "optim")))$tTable
    Chat[i,] = tSummary[1:6, 1]
    SEmtx[i,] = tSummary[1:6, 2]
    DFmtx[i,] = tSummary[1:6, 3]
  }, error=function(e){})
}

colnames(Chat) = colnames(SEmtx) = colnames(DFmtx) = c('AAD', 'ASD', 'BD', 'CTL','MDD','SCZ')
saveRDS(list(Chat = Chat, SE = SEmtx, DF = DFmtx), 'data/results/Microarray_compiledGLM.rds')
