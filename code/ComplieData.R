asd_meta = readRDS('data/results/control/Microarray_ASD_metaanalysis.rds')
scz_meta = readRDS('data/results/control/Microarray_SCZ_metaanalysis.rds')
bd_meta = readRDS('data/results/control/Microarray_BD_metaanalysis.rds')
mdd_meta = readRDS('data/results/control/Microarray_MDD_metaanalysis.rds')
aad_meta = readRDS('data/results/control/Microarray_AAD_metaanalysis.rds')

multi= vector(mode="list", length=5)
multi[[1]]$beta = asd_meta$beta;  multi[[1]]$se = asd_meta$SE;
multi[[1]]$DF = asd_meta$DF;  multi[[1]]$p.value = asd_meta$p;
multi[[1]]$genes = row.names(asd_meta)

multi[[2]]$beta = scz_meta$beta;  multi[[2]]$se = scz_meta$SE;
multi[[2]]$DF = scz_meta$DF;  multi[[2]]$p.value = scz_meta$p;
multi[[2]]$genes = row.names(scz_meta)

multi[[3]]$beta = bd_meta$beta;  multi[[3]]$se = bd_meta$SE;
multi[[3]]$DF = bd_meta$DF;  multi[[3]]$p.value = bd_meta$p;
multi[[3]]$genes = row.names(bd_meta)

multi[[4]]$beta = mdd_meta$beta;  multi[[4]]$se = mdd_meta$SE;
multi[[4]]$DF = mdd_meta$DF;  multi[[4]]$p.value = mdd_meta$p;
multi[[4]]$genes = row.names(mdd_meta)

multi[[5]]$beta = aad_meta$beta;  multi[[5]]$se = aad_meta$SE;
multi[[5]]$DF = aad_meta$DF;  multi[[5]]$p.value = aad_meta$p;
multi[[5]]$genes = row.names(aad_meta)

names(multi) = c("ASD", "SCZ", "BD", "MDD", "AAD")

genes = intersect(multi[[1]]$genes, multi[[2]]$genes)
for(i in 3:5){genes = intersect(genes, multi[[i]]$genes)}

for(i in 1:5){
  ind = match(genes, multi[[i]]$genes)
  multi[[i]]$beta = multi[[i]]$beta[ind]
  multi[[i]]$se = multi[[i]]$se[ind]
  multi[[i]]$DF = multi[[i]]$DF[ind]
  multi[[i]]$p.value = multi[[i]]$p.value[ind]
  multi[[i]]$genes = multi[[i]]$genes[ind]
}

beta = cbind(multi[[1]]$beta, multi[[2]]$beta)
se = cbind(multi[[1]]$se, multi[[2]]$se)
df = cbind(multi[[1]]$DF, multi[[2]]$DF)
p = cbind(multi[[1]]$p.value, multi[[2]]$p.value)
for(i in 3:5){
  beta = cbind(beta, multi[[i]]$beta)
  se = cbind(se, multi[[i]]$se)
  df = cbind(df, multi[[i]]$DF)
  p = cbind(p, multi[[i]]$p.value)
}
colnames(beta) = colnames(se) = colnames(df) = colnames(p) = c("ASD", "SCZ", "BD", "MDD", "AAD")
row.names(beta) = row.names(se) = row.names(df) = row.names(p) = genes
saveRDS(list(beta = beta, se=se, df=df, p = p), 'data/results/control/CompiledData.rds')
