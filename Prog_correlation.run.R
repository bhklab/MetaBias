
dir.create("/results/Cor")
####################################################################
## breast data: correlation expression data
####################################################################

load("/data/breast/DrugMeta.breast.zfixsample.RData")
load("/data/breast/DrugMeta.breast.fixsample.RData")

id <- z.exprs.ccle.data$gene_name
z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data, z.exprs.gray.data, 
                z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data, z.exprs.uhn.data)

study <- c("ccle", "ctrpv", "gcsi", "gray", "gdsc1", "gdsc2", "uhn")


cor.genes <- lapply(1:length(id), function(k){
  
cor.dat <- lapply(1:length(z.exprs), function(i){
    
   dat <- z.exprs[[i]]    
   as.numeric(dat[dat$gene_name == id[k], -1])
     
  })
  
  res <- do.call(cbind, cor.dat)
  cor.res <- cor(res, use = "complete.obs", method = "pearson")
  cor.val <- cor.res[upper.tri(cor.res, diag = FALSE)]
  mean.cor <- mean(cor.val)
  median.cor <- median(cor.val)
  data.frame(gene = id[k], mean.cor, median.cor)
  
})

cor.res <- do.call(rbind, cor.genes)

save(cor.res, file="/results/Cor/cor.breast.RData")

##################################################################################
## breast data: correlation expression data
##################################################################################

load("/data/pan/DrugMeta.pan.zfixsample.RData")
load("/data/pan/DrugMeta.pan.fixsample.mod.RData")

id <- z.exprs.ccle.data$gene_name
z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data,
                z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data)

study <- c("ccle", "ctrpv", "gcsi", "gdsc1", "gdsc2")


cor.genes <- lapply(1:length(id), function(k){
  
  cor.dat <- lapply(1:length(z.exprs), function(i){
    
    dat <- z.exprs[[i]]    
    as.numeric(dat[dat$gene_name == id[k], -1])
    
  })
  
  res <- do.call(cbind, cor.dat)
  cor.res <- cor(res, use = "complete.obs", method = "pearson")
  cor.val <- cor.res[upper.tri(cor.res, diag = FALSE)]
  mean.cor <- mean(cor.val)
  median.cor <- median(cor.val)
  data.frame(gene = id[k], mean.cor, median.cor)
  
})

cor.res <- do.call(rbind, cor.genes)

save(cor.res, file="/results/Cor/cor.pan.RData")


