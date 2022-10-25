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
method.meta <- c("Jeffreys")

combine.lm.res <- lapply(1:length(drugs), function(k){
  
  res <- lapply(1:length(lm.dat), function(i){
    
     data.frame(study.id = study[i], lm.dat[[i]][lm.dat[[i]]$drug == drugs[k], ]) 

      })
  
  do.call(rbind, res)
  
})

combine.lm.res <- do.call(rbind, combine.lm.res)
gene.id <- unique(combine.lm.res$gene.name)

###################################################################
## combine effect sizes
###################################################################
meta.effect.size <- lapply(1:length(drugs), function(k){ 
  
  print(k)
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
    res.meta.gene <- lapply(1:length(method.meta), function(j){ 
      
      res.meta <- lapply(1:length(gene.id), function(i){ # length(gene.id)
      
      sub.dat <- sub.res[sub.res$gene.name == gene.id[i], ]
      meta.bayes.fun(sub.dat$estimate, sub.dat$std, method.meta[j], sub.dat$study.id)
      
    })
    
    res.meta.method <- do.call(rbind, res.meta)
    res.meta.method$method.meta <- method.meta[j]
    #res.meta.method$padj <- p.adjust(res.meta.method$pval.random, method = "BH")
    data.frame(gene.name = gene.id,
               drug = drugs[k],
               res.meta.method)

  })
  
  res.meta.gene.drug <-  do.call(rbind, res.meta.gene)

  })

meta.effect.size.res <- do.call(rbind, meta.effect.size)

save(meta.effect.size.res, 
     file = "/results/Meta/Breast/meta.breast.bayes.RData")

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
method.meta <- c("Jeffreys")

combine.lm.res <- lapply(1:length(drugs), function(k){
  
  res <- lapply(1:length(lm.dat), function(i){
    
    data.frame(study.id = study[i], lm.dat[[i]][lm.dat[[i]]$drug == drugs[k], ]) 
    
  })
  
  do.call(rbind, res)
  
})

combine.lm.res <- do.call(rbind, combine.lm.res)
gene.id <- unique(combine.lm.res$gene.name)

###################################################################
## combine effect sizes
###################################################################
meta.effect.size <- lapply(1:length(drugs), function(k){
  
  print(k)
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res.meta.gene <- lapply(1:length(method.meta), function(j){
    
    res.meta <- lapply(1:length(gene.id), function(i){ #length(gene.id)
      
      sub.dat <- sub.res[sub.res$gene.name == gene.id[i], ]
      meta.bayes.fun(sub.dat$estimate, sub.dat$std, method.meta[j], sub.dat$study.id)
      
    })
    
    res.meta.method <- do.call(rbind, res.meta)
    res.meta.method$method.meta <- method.meta[j]
    #res.meta.method$padj <- p.adjust(res.meta.method$pval.random, method = "BH")
    data.frame(gene.name = gene.id,
               drug = drugs[k],
               res.meta.method)
    
  })
  
  res.meta.gene.drug <-  do.call(rbind, res.meta.gene)
  
})

meta.effect.size.res <- do.call(rbind, meta.effect.size)

save(meta.effect.size.res, 
     file = "/results/Meta/Pan/meta.pan.bayes.RData")


