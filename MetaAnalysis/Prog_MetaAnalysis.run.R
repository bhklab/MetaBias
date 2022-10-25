## Packages

source("/code/MetaAnalysis/Prog_metaAnalysis_fun.R")
source("/code/LinearModel/Prog_LinearModel.run.R")

dir.create("/results/Meta")
dir.create("/results/Meta/Breast")
dir.create("/results/Meta/Pan")
#############################################################
#############################################################
## Load Breast linear model results 
#############################################################
#############################################################

load("/results/LM/Breast/lm.breast.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Breast/lm.breast.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Breast/lm.breast.gray.RData")
res.gray <- lm.res
load("/results/LM/Breast/lm.breast.uhn.RData")
res.uhn <- lm.res
load("/results/LM/Breast/lm.breast.ctrpv.RData")
res.ctrpv <- lm.res
load("/results/LM/Breast/lm.breast.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Breast/lm.breast.gdsc2.RData")
res.gdsc2 <- lm.res

###########################################################
## merge data
###########################################################

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel") 
lm.dat <- list(res.ccle, res.ctrpv, res.gcsi, res.gray, res.gdsc1, res.gdsc2, res.uhn)
study <- c("ccle", "ctrp", "gcsi", "gray", "gdsc1", "gdsc2", "uhn")
method.meta <- c("DL", "SJ", "HS", "EB", "HE", "PM")

combine.lm.res <- lapply(1:length(drugs), function(k){
  
  res <- lapply(1:length(lm.dat), function(i){
    
     data.frame(study.id = study[i], lm.dat[[i]][lm.dat[[i]]$drug == drugs[k], ]) 

      })
  
  do.call(rbind, res)
  
})

combine.lm.res <- do.call(rbind, combine.lm.res)

gene.id <- unique(combine.lm.res$gene.name)

####################################################################
## combine p values
####################################################################
meta.pvalue <- lapply(1:length(drugs), function(k){
  
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res <- lapply(1:length(gene.id), function(i){
    
      fisher.p <- sumLog(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"])
      stouffer.p <- sumZ(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"]) 
      acat.p <- ACAT(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"], 
                       weights = NULL, is.check = TRUE)
      hmp.p <- hmp.stat(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"])
      
      data.frame(gene.name = gene.id[i],
                        drug = drugs[k],
                        fisher.p = fisher.p$p, 
                        stouffer.p = stouffer.p$p,
                        acat.p = acat.p,
                        hmp.p  = hmp.p)
  })
  
  meta.res <- do.call(rbind, res)
  meta.res$padj.fisher <- p.adjust(meta.res$fisher.p, method = "BH")
  meta.res$padj.stouffer <- p.adjust(meta.res$stouffer.p, method = "BH")
  meta.res$padj.acat <- p.adjust(meta.res$acat.p, method = "BH")
  meta.res$padj.hmp <- p.adjust(meta.res$hmp.p, method = "BH")
  meta.res
  
})

meta.pvalue.res <- do.call(rbind, meta.pvalue)
rownames(meta.pvalue.res) <- NULL
###################################################################
## combine effect sizes
###################################################################
meta.effect.size <- lapply(1:length(drugs), function(k){ 
  
  print(k)
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
    res.meta.gene <- lapply(1:length(method.meta), function(j){ 
      
      res.meta <- lapply(1:length(gene.id), function(i){
      
      sub.dat <- sub.res[sub.res$gene.name == gene.id[i], ]
      meta.fun(sub.dat$estimate, sub.dat$std, method = method.meta[j], sub.dat$study.id)
      
    })
    
    res.meta.method <- do.call(rbind, res.meta)
    res.meta.method$method.meta <- method.meta[j]
    res.meta.method$padj <- p.adjust(res.meta.method$pval.random, method = "BH")
    data.frame(gene.name = gene.id,
               drug = drugs[k],
               res.meta.method)

  })
  
  res.meta.gene.drug <-  do.call(rbind, res.meta.gene)

  })

meta.effect.size.res <- do.call(rbind, meta.effect.size)

save(meta.effect.size.res, meta.pvalue.res, 
     file = "/results/Meta/Breast/meta.breast.RData")


#############################################################
#############################################################
## Load Pan linear model results 
#############################################################
#############################################################

load("/results/LM/Pan/lm.pan.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Pan/lm.pan.ctrpv.RData")
res.ctrpv <- lm.res
load("/results/LM/Pan/lm.pan.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Pan/lm.pan.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Pan/lm.pan.gdsc2.RData")
res.gdsc2 <- lm.res
###########################################################
## merge data
###########################################################

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel") 
lm.dat <- list(res.ccle, res.ctrpv, res.gcsi, res.gdsc1, res.gdsc2)
study <- c("ccle", "ctrp", "gcsi", "gdsc1", "gdsc2")
method.meta <- c("DL", "SJ", "HS", "EB", "HE", "PM")

combine.lm.res <- lapply(1:length(drugs), function(k){
  
  res <- lapply(1:length(lm.dat), function(i){
    
    data.frame(study.id = study[i], lm.dat[[i]][lm.dat[[i]]$drug == drugs[k], ]) 
    
  })
  
  do.call(rbind, res)
  
})

combine.lm.res <- do.call(rbind, combine.lm.res)

gene.id <- unique(combine.lm.res$gene.name)
####################################################################
## combine p values
####################################################################
meta.pvalue <- lapply(1:length(drugs), function(k){
  
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res <- lapply(1:length(gene.id), function(i){
    
    fisher.p <- sumLog(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"])
    stouffer.p <- sumZ(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"]) 
    acat.p <- ACAT(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"], 
                   weights = NULL, is.check = TRUE)
    hmp.p <- hmp.stat(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"])
    data.frame(gene.name = gene.id[i],
               drug = drugs[k],
               fisher.p = fisher.p$p, 
               stouffer.p = stouffer.p$p,
               acat.p = acat.p,
               hmp.p  = hmp.p)
  })
  
  meta.res <- do.call(rbind, res)
  meta.res$padj.fisher <- p.adjust(meta.res$fisher.p, method = "BH")
  meta.res$padj.stouffer <- p.adjust(meta.res$stouffer.p, method = "BH")
  meta.res$padj.acat <- p.adjust(meta.res$acat.p, method = "BH")
  meta.res$padj.hmp <- p.adjust(meta.res$hmp.p, method = "BH")
  meta.res
  
})

meta.pvalue.res <- do.call(rbind, meta.pvalue)
rownames(meta.pvalue.res) <- NULL
###################################################################
## combine effect sizes
###################################################################
meta.effect.size <- lapply(1:length(drugs), function(k){
  
  print(k)
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res.meta.gene <- lapply(1:length(method.meta), function(j){
    
    res.meta <- lapply(1:length(gene.id), function(i){
      
      sub.dat <- sub.res[sub.res$gene.name == gene.id[i], ]
      meta.fun(sub.dat$estimate, sub.dat$std, method.meta[j], sub.dat$study.id)
      
    })
    
    res.meta.method <- do.call(rbind, res.meta)
    res.meta.method$method.meta <- method.meta[j]
    res.meta.method$padj <- p.adjust(res.meta.method$pval.random, method = "BH")
    data.frame(gene.name = gene.id,
               drug = drugs[k],
               res.meta.method)
    
  })
  
  res.meta.gene.drug <-  do.call(rbind, res.meta.gene)
  
})

meta.effect.size.res <- do.call(rbind, meta.effect.size)

save(meta.effect.size.res, meta.pvalue.res, 
     file = "/results/Meta/Pan/meta.pan.RData")


