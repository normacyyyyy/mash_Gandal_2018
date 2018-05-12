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

asd_meta = matrix(NA, nrow=nrow(ASD_datExpr), ncol=4)
for(i in 1:nrow(ASD_datExpr)) {
  if(i%%100==0) print(i)
  expr = ASD_datExpr[i,]
  tryCatch({
    asd_meta[i,] = summary(lme(expr~ Group + Study,data = ASD_datMeta, random=~1|Subject, control = list(opt = "optim")))$tTable[2,c(1,2,3,5)]
  }, error=function(e){})
}
asd_meta = as.data.frame(asd_meta)
colnames(asd_meta) = c("beta", "SE", "DF", "p")
rownames(asd_meta) = genes
asd_meta$fdr = p.adjust(asd_meta$p, "fdr")
asd_meta$symbol = ASD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, ASD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
to_keep = !apply(is.na(asd_meta),1,any)
asd_meta = asd_meta[to_keep,]
saveRDS(asd_meta, 'data/results/control/Microarray_ASD_metaanalysis.rds')


# SCZ
files = dir("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("SCZ",files)]
n = length(files)
SCZ_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  SCZ_multiExpr[[i]]$datExpr = datExpr
  SCZ_multiExpr[[i]]$datMeta = datMeta
  SCZ_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(rownames(SCZ_multiExpr[[1]]$datExpr), rownames(SCZ_multiExpr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(SCZ_multiExpr[[i]]$datExpr))
for(i in 1:n) SCZ_multiExpr[[i]]$datExpr = SCZ_multiExpr[[i]]$datExpr[match(genes,rownames(SCZ_multiExpr[[i]]$datExpr)),]

SCZ_datExpr = data.frame(row.names = genes)
SCZ_datMeta=data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  SCZ_datExpr = cbind(SCZ_datExpr, SCZ_multiExpr[[i]]$datExpr)

  if("Group.SCZ" %in% colnames(SCZ_multiExpr[[i]]$datMeta)) SCZ_multiExpr[[i]]$datMeta$Group = SCZ_multiExpr[[i]]$datMeta$Group.SCZ   #Only for unique controls
  SCZ_datMeta = rbind(SCZ_datMeta, SCZ_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group")])

}
SCZ_datMeta = SCZ_datMeta[-1,]
to_keep = SCZ_datMeta$Group %in% c("CTL", "SCZ")
SCZ_datExpr = SCZ_datExpr[,to_keep]; SCZ_datMeta = SCZ_datMeta[to_keep,]
SCZ_datMeta$Group = factor(SCZ_datMeta$Group)
SCZ_datMeta$Study = as.factor(SCZ_datMeta$Study)
SCZ_datMeta$Subject = as.factor(SCZ_datMeta$Subject)

scz_meta = matrix(NA, nrow=nrow(SCZ_datExpr), ncol=4)
for(i in 1:nrow(SCZ_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(SCZ_datExpr[i,])
  tryCatch({
    scz_meta[i,] = summary(lme(expr~ Group + Study,data = SCZ_datMeta, random=~1|Subject))$tTable[2,c(1,2,3,5)]
  }, error=function(e){})
}

scz_meta=as.data.frame(scz_meta)
colnames(scz_meta) = c("beta", "SE", "DF", "p")
rownames(scz_meta) = genes
scz_meta$fdr = p.adjust(scz_meta$p, "fdr")
scz_meta$symbol=SCZ_multiExpr[[1]]$datProbes$external_gene_id[match(genes, SCZ_multiExpr[[1]]$datProbes$ensembl_gene_id)]
scz_meta= scz_meta[!apply(is.na(scz_meta),1,any),]
saveRDS(scz_meta, 'data/results/control/Microarray_SCZ_metaanalysis.rds')

#BD
files = dir("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("_BD_",files)]
n=length(files)
BD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  BD_multiExpr[[i]]$datExpr = datExpr
  BD_multiExpr[[i]]$datMeta = datMeta
  BD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(rownames(BD_multiExpr[[1]]$datExpr), rownames(BD_multiExpr[[2]]$datExpr))
for(i in 3:n) genes = intersect(genes, rownames(BD_multiExpr[[i]]$datExpr))
for(i in 1:n) BD_multiExpr[[i]]$datExpr = BD_multiExpr[[i]]$datExpr[match(genes,rownames(BD_multiExpr[[i]]$datExpr)),]

BD_datExpr = data.frame(row.names = genes)
BD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  BD_datExpr = cbind(BD_datExpr, BD_multiExpr[[i]]$datExpr)
  if("Group.BD" %in% colnames(BD_multiExpr[[i]]$datMeta)) BD_multiExpr[[i]]$datMeta$Group = BD_multiExpr[[i]]$datMeta$Group.BD ## --> only for unique controls
  BD_datMeta = rbind(BD_datMeta, BD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group")])
}
BD_datMeta = BD_datMeta[-1,]
to_keep = BD_datMeta$Group %in% c("CTL", "BD")
BD_datExpr = BD_datExpr[,to_keep]; BD_datMeta = BD_datMeta[to_keep,]
BD_datMeta$Group = factor(BD_datMeta$Group, levels=c("CTL", "BD"))
BD_datMeta$Study = as.factor(BD_datMeta$Study)
BD_datMeta$Subject = as.factor(BD_datMeta$Subject)

bd_meta = matrix(NA, nrow=nrow(BD_datExpr), ncol=4)
for(i in 1:nrow(BD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(BD_datExpr[i,])
  tryCatch({
    bd_meta[i,] = summary(lme(expr~ Group + Study,data = BD_datMeta, random=~1|Subject))$tTable[2,c(1,2,3,5)]
  }, error=function(e){})
}

bd_meta=as.data.frame(bd_meta)
colnames(bd_meta) = c("beta", "SE", "DF", "p")
rownames(bd_meta) = genes
bd_meta$fdr = p.adjust(bd_meta$p, "fdr")
bd_meta$symbol=BD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, BD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
bd_meta= bd_meta[!apply(is.na(bd_meta),1,any),]

saveRDS(bd_meta, "data/results/control/Microarray_BD_metaanalysis.rds")

# MDD
files = dir("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", pattern="_CR_"); files = files[grep("_MDD_",files)]
n=length(files)
MDD_multiExpr = vector(mode="list",length = n)
for(i in 1:n) {
  load(paste("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/", files[[i]], sep=""))
  MDD_multiExpr[[i]]$datExpr = datExpr
  MDD_multiExpr[[i]]$datMeta = datMeta
  MDD_multiExpr[[i]]$datProbes = datProbes
  rm(datExpr,datMeta,datProbes)
}

genes= intersect(rownames(MDD_multiExpr[[1]]$datExpr), rownames(MDD_multiExpr[[2]]$datExpr))
for(i in 1:n) MDD_multiExpr[[i]]$datExpr = MDD_multiExpr[[i]]$datExpr[match(genes,rownames(MDD_multiExpr[[i]]$datExpr)),]

MDD_datExpr = data.frame(row.names = genes)
MDD_datMeta=data.frame(Study=NA, Subject=NA, Group=NA)
for(i in 1:n) {
  MDD_datExpr = cbind(MDD_datExpr, MDD_multiExpr[[i]]$datExpr)
  MDD_datMeta = rbind(MDD_datMeta, MDD_multiExpr[[i]]$datMeta[,c("Study", "Subject","Group")])
}
MDD_datMeta = MDD_datMeta[-1,]
to_keep = MDD_datMeta$Group %in% c("CTL", "MDD")
MDD_datExpr = MDD_datExpr[,to_keep]; MDD_datMeta = MDD_datMeta[to_keep,]
MDD_datMeta$Group = factor(MDD_datMeta$Group, levels=c("CTL", "MDD"))
MDD_datMeta$Study = as.factor(MDD_datMeta$Study)
MDD_datMeta$Subject = as.factor(MDD_datMeta$Subject)

mdd_meta = matrix(NA, nrow=nrow(MDD_datExpr), ncol=4)
for(i in 1:nrow(MDD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(MDD_datExpr[i,])
  tryCatch({
    mdd_meta[i,] = summary(lme(expr~ Group + Study,data = MDD_datMeta, random=~1|Subject))$tTable[2,c(1,2,3,5)]
  }, error=function(e){})
}

mdd_meta=as.data.frame(mdd_meta)
colnames(mdd_meta) = c("beta", "SE", "DF", "p")
rownames(mdd_meta) = genes
mdd_meta$fdr = p.adjust(mdd_meta$p, "fdr")
mdd_meta$symbol=MDD_multiExpr[[1]]$datProbes$external_gene_id[match(genes, MDD_multiExpr[[1]]$datProbes$ensembl_gene_id)]
mdd_meta= mdd_meta[!apply(is.na(mdd_meta),1,any),]
saveRDS(mdd_meta, 'data/results/control/Microarray_MDD_metaanalysis.rds')

##AAD
load("data/working_data/Microarray/03_NormalizedBalanced_ComBat_CR_cleaned/Microarray_AAD_mayfield_normalized_CR_regressed.RData")
AAD_datExpr = datExpr; rm(datExpr)
AAD_datMeta = datMeta; rm(datMeta)
AAD_datProbes = datProbes; rm(datProbes)
genes= rownames(AAD_datExpr)
aad_meta = matrix(NA, nrow=nrow(AAD_datExpr), ncol=4)
for(i in 1:nrow(AAD_datExpr)) {
  if(i%%100==0) print(i)
  expr = as.numeric(AAD_datExpr[i,])
  tryCatch({
    aad_meta[i,] = summary(lme(expr~ Group,data = AAD_datMeta, random=~1|Subject))$tTable[2,c(1,2,3,5)]
  }, error=function(e){})
}
aad_meta=as.data.frame(aad_meta)
colnames(aad_meta) = c("beta", "SE", "DF", "p")
rownames(aad_meta) = genes
aad_meta$fdr = p.adjust(aad_meta$p, "fdr")
aad_meta$symbol=AAD_datProbes$external_gene_id[match(genes, AAD_datProbes$ensembl_gene_id)]
aad_meta= aad_meta[!apply(is.na(aad_meta),1,any),]
saveRDS(aad_meta, 'data/results/control/Microarray_AAD_metaanalysis.rds')
