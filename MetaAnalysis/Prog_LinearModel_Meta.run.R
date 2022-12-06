## Packages

source("/code/Prog_LinearModel_Meta_fun.R")
source("/code/Prog_upset_meta.R")

dir.create("/results/LM")
dir.create("/results/LM/Breast")
dir.create("/results/LM/Pan")

dir.create("/results/Meta")
dir.create("/results/Meta/Breast")
dir.create("/results/Meta/Pan")

dir.create("/results/MetaVisualization")
dir.create("/results/MetaVisualization/Breast")
dir.create("/results/MetaVisualization/Pan")

################################################################################################
################################################################################################
############################### (adjusted) linear model fitting ################################
################################################################################################
################################################################################################

## linear model: breast cancer cell line

load("/data/breast/DrugMeta.breast.zfixsample.RData")
load("/data/breast/DrugMeta.breast.fixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

z.aac <- list(mod.aac.ccle, mod.aac.ctrpv, mod.aac.gcsi,  mod.aac.gray, mod.aac.gdsc.v1, mod.aac.gdsc.v2, mod.aac.uhn)

z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data, z.exprs.gray.data, z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data, z.exprs.uhn.data)

study <- c("ccle", "ctrpv", "gcsi", "gray", "gdsc1", "gdsc2", "uhn")

for(i in 1:length(study)){
  
  print(i)
  lm.res <- lm.fun(exprs = z.exprs[[i]], aac= z.aac[[i]], drug = drugs)   
  save(lm.res, file = paste("/results/LM/Breast", 
                                     paste("lm.breast", study[i], "RData", sep="."), sep="/") ) 
  
}


## linear model: pan cancer cell line

load("/data/pan/DrugMeta.pan.fixsample.mod.RData")
load("/data/pan/DrugMeta.pan.zfixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

z.aac <- list(mod.aac.ccle, mod.aac.ctrpv, mod.aac.gcsi,
              mod.aac.gdsc.v1, mod.aac.gdsc.v2)

z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data, 
                z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data)

ccl.dat <- list(pan.ccl.ccle, pan.ccl.ctrpv, pan.ccl.gcsi, pan.ccl.gdsc.v1, pan.ccl.gdsc.v2)

study <- c("ccle", "ctrpv", "gcsi", "gdsc1", "gdsc2")

for(i in 1:length(study)){
  
  print(i)
  lm.res <- lm.adjusted.fun(exprs = z.exprs[[i]], aac= z.aac[[i]], drug = drugs, ccl = ccl.dat[[i]])   
  save(lm.res, file = paste("/results/LM/Pan", 
                            paste("lm.pan", study[i], "RData", sep="."), sep="/") ) 
  
}

################################################################################################
################################################################################################
############################### meta-analysis fitting ##########################################
################################################################################################
################################################################################################

## Load Breast linear model results 

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

## merge linear model findings across studies

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

## combine p values

meta.pvalue <- lapply(1:length(drugs), function(k){
  
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res <- lapply(1:length(gene.id), function(i){
    
    fisher.p <- sumLog(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"])
    stouffer.p <- sumZ(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"]) 
    
    data.frame(gene.name = gene.id[i],
               drug = drugs[k],
               fisher.p = fisher.p$p, 
               stouffer.p = stouffer.p$p)
  })
  
  meta.res <- do.call(rbind, res)
  meta.res$padj.fisher <- p.adjust(meta.res$fisher.p, method = "BH")
  meta.res$padj.stouffer <- p.adjust(meta.res$stouffer.p, method = "BH")
  
  meta.res
  
})

meta.pvalue.res <- do.call(rbind, meta.pvalue)
rownames(meta.pvalue.res) <- NULL


## combine effect sizes

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



## Load Pan linear model results 

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


## merge linear model findings across studies

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

## combine p values

meta.pvalue <- lapply(1:length(drugs), function(k){
  
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res <- lapply(1:length(gene.id), function(i){
    
    fisher.p <- sumLog(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"])
    stouffer.p <- sumZ(sub.res[sub.res$gene.name == gene.id[i], "Pvalue"]) 
    
    data.frame(gene.name = gene.id[i],
               drug = drugs[k],
               fisher.p = fisher.p$p, 
               stouffer.p = stouffer.p$p)
  })
  
  meta.res <- do.call(rbind, res)
  meta.res$padj.fisher <- p.adjust(meta.res$fisher.p, method = "BH")
  meta.res$padj.stouffer <- p.adjust(meta.res$stouffer.p, method = "BH")
 
  meta.res
  
})

meta.pvalue.res <- do.call(rbind, meta.pvalue)
rownames(meta.pvalue.res) <- NULL

## combine effect sizes

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

####################################################################################################
####################################################################################################
#################################### meta-analysis: Bayesian approach ##############################
####################################################################################################
####################################################################################################

## Load Breast linear model results 

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

## merge linear model findings across studies

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

## combine effect sizes: Bayesian approach

# NOTE: change the loop to length(gene.id)
set.seed(135)
gene.id <- c(sample(gene.id, 50),
             "CHMP7", "FBXO3", "DUOX2", "EGFR", "CD63", "ZNF143",
             "PSMB3", "C3orf52", "PITX3", "ERBB2", "CD63", "CD84",
             "IGF2BP3", "TGM3", "ABCG4", "S100A1", "COL11A1", "SCRT1")

meta.effect.size <- lapply(1:length(drugs), function(k){ 
  
  print(k)
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res.meta.gene <- lapply(1:length(method.meta), function(j){ 
    
    res.meta <- lapply(1:length(gene.id), function(i){ 
      
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


## Load Pan linear model results 

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

## merge linear model findings across studies

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


## combine effect sizes

# NOTE: change the loop to length(gene.id)
set.seed(135)
gene.id <- c(sample(gene.id, 50),
             "C1orf116", "SPRR1B", "KRT1", "EGFR", "COX7A2", "ROS1" ,
             "TDRD1", "ZNF365", "KRT1", "ERBB2", "MAPK7", "CYP2A13", 
             "C4orf19", "ZNF688", "CPSF6", "S100A1", "NOC3L", "CD244")

meta.effect.size <- lapply(1:length(drugs), function(k){
  
  print(k)
  sub.res <- combine.lm.res[combine.lm.res$drug == drugs[k], ]
  
  res.meta.gene <- lapply(1:length(method.meta), function(j){
    
    res.meta <- lapply(1:length(gene.id), function(i){
      
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

##################################################################################################
##################################################################################################
###################################### visualization plot ########################################
##################################################################################################
##################################################################################################

####################################################################
## UpSet plot to compare different meta-analysis methods: breast data
####################################################################

load("/results/Meta/Breast/meta.breast.RData")

id <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(i in 1:length(id)){
  
  res.adj.Fisher <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.fisher")]
  res.adj.Stouffer <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.stouffer")]
 
  res.adj.DL <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "DL", c("gene.name","padj") ]
  res.adj.EB <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "EB", c("gene.name","padj") ]
  res.adj.HE <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "HE", c("gene.name","padj") ]
  res.adj.HS <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "HS", c("gene.name","padj") ]
  res.adj.PM <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "PM", c("gene.name","padj") ]
  res.adj.SJ <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "SJ", c("gene.name","padj") ]
  
  sig.fisher <- res.adj.Fisher[res.adj.Fisher$padj.fisher < 0.05, "gene.name"]
  sig.stouffer <- res.adj.Stouffer[res.adj.Stouffer$padj.stouffer < 0.05, "gene.name"]
  
  sig.DL <- res.adj.DL[res.adj.DL$padj < 0.05, "gene.name"]
  sig.SJ <- res.adj.SJ[res.adj.SJ$padj < 0.05, "gene.name"]
  sig.HS <- res.adj.HS[res.adj.HS$padj < 0.05, "gene.name"]
  sig.HE <- res.adj.HE[res.adj.HE$padj < 0.05, "gene.name"]
  sig.EB <- res.adj.EB[res.adj.EB$padj < 0.05, "gene.name"]
  sig.PM <- res.adj.PM[res.adj.PM$padj < 0.05, "gene.name"]
  
  listInput <- list(DL = sig.DL,
                    SJ = sig.SJ,
                    HS = sig.HS, 
                    HE = sig.HE,
                    EB = sig.EB,
                    PM = sig.PM,
                    Fisher = sig.fisher,
                    Stouffer = sig.stouffer)
  
  
  dir0 <- "/results/Meta/Breast"
  dir <- paste(paste(dir0, paste("SFig3_cont", id[i], sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    20, #length
    15, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  col.id <- c("maroon","blue2", "darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3", "salmon3")
  sets = c("DL", "SJ", "HS", "HE", "EB", "PM", "Fisher", "Stouffer")
  
  sets.bar.color= col.id
  
  p <- upset.fun(fromList(listInput), paste(id[i], "FDR < 0.05", sep=", "))
  print(p)
  
  dev.off()
  
}

######################################################################
## UpSet plot to compare different meta-analysis methods: pan data
######################################################################

load("/results/Meta/Pan/meta.pan.RData")

id <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(i in 1:length(id)){
  
  res.adj.Fisher <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.fisher")]
  res.adj.Stouffer <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.stouffer")]
  
  res.adj.DL <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "DL", c("gene.name","padj") ]
  res.adj.EB <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "EB", c("gene.name","padj") ]
  res.adj.HE <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "HE", c("gene.name","padj") ]
  res.adj.HS <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "HS", c("gene.name","padj") ]
  res.adj.PM <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "PM", c("gene.name","padj") ]
  res.adj.SJ <- meta.effect.size.res[meta.effect.size.res$drug == id[i] &
                                       meta.effect.size.res$method.meta == "SJ", c("gene.name","padj") ]
  
  sig.fisher <- res.adj.Fisher[res.adj.Fisher$padj.fisher < 0.05, "gene.name"]
  sig.stouffer <- res.adj.Stouffer[res.adj.Stouffer$padj.stouffer < 0.05, "gene.name"]
 
  sig.DL <- res.adj.DL[res.adj.DL$padj < 0.05, "gene.name"]
  sig.SJ <- res.adj.SJ[res.adj.SJ$padj < 0.05, "gene.name"]
  sig.HS <- res.adj.HS[res.adj.HS$padj < 0.05, "gene.name"]
  sig.HE <- res.adj.HE[res.adj.HE$padj < 0.05, "gene.name"]
  sig.EB <- res.adj.EB[res.adj.EB$padj < 0.05, "gene.name"]
  sig.PM <- res.adj.PM[res.adj.PM$padj < 0.05, "gene.name"]
  
  listInput <- list(DL = sig.DL,
                    SJ = sig.SJ,
                    HS = sig.HS, 
                    HE = sig.HE,
                    EB = sig.EB,
                    PM = sig.PM,
                    Fisher = sig.fisher,
                    Stouffer = sig.stouffer)
  
  
  dir0 <- "/results/Meta/Pan"
  dir <- paste(paste(dir0, paste("SFig3_cont", id[i], sep="_"), sep="/"), "jpeg", sep=".")
  Cairo::Cairo(
    20, #length
    15, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  col.id <- c("maroon","blue2", "darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3", "salmon3")
  sets = c("DL", "SJ", "HS", "HE", "EB", "PM", "Fisher", "Stouffer")
  
  sets.bar.color= col.id
  
  p <- upset.fun(fromList(listInput), paste(id[i], "FDR < 0.05", sep=", "))
  print(p)
  
  dev.off()
  
}

####################################################################
## forest plot: breast data
####################################################################

load("/results/LM/Breast/lm.breast.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Breast/lm.breast.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Breast/lm.breast.gray.RData")
res.gray <- lm.res
load("/results/LM/Breast/lm.breast.uhn.RData")
res.uhn <- lm.res
load("/results/LM/Breast/lm.breast.ctrpv.RData")
res.ctrp <- lm.res
load("/results/LM/Breast/lm.breast.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Breast/lm.breast.gdsc2.RData")
res.gdsc2 <- lm.res

int <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id <- c("EGFR", "ERBB2", "S100A1")

for(k in 1:length(int)){
  
  sub.res.ccle <- res.ccle[res.ccle$gene.name == id[k] & 
                             res.ccle$drug == int[k], ]
  sub.res.gray <- res.gray[res.gray$gene.name == id[k] & 
                             res.gray$drug == int[k], ]
  sub.res.uhn <- res.uhn[res.uhn$gene.name == id[k] & 
                           res.uhn$drug == int[k], ]
  sub.res.gcsi <- res.gcsi[res.gcsi$gene.name == id[k] & 
                             res.gcsi$drug == int[k], ]
  sub.res.gdsc1 <- res.gdsc1[res.gdsc1$gene.name == id[k] & 
                               res.gdsc1$drug == int[k], ]
  sub.res.gdsc2 <- res.gdsc2[res.gdsc2$gene.name == id[k] & 
                               res.gdsc2$drug == int[k], ]
  sub.res.ctrp <- res.ctrp[res.ctrp$gene.name == id[k] & 
                             res.ctrp$drug == int[k], ]
  
  res <- rbind(sub.res.ccle, sub.res.ctrp, sub.res.gcsi,
               sub.res.gdsc1, sub.res.gdsc2, sub.res.gray, sub.res.uhn)
  study.id <- c("CCLE", "CTRP", "gCSI", "GDSC1","GDSC2", "GRAY", "UHNBreast")
  df.id <- paste("n", res$df + 2, sep=" = ")
  rownames(res) <- paste(study.id, df.id, sep=", ")
  
  x <-  res$estimate
  y <-  res$std
  dat.xy <- data.frame(x, y)
  rownames(dat.xy) <- rownames(res)
  res <- metagen(x, y, studlab = rownames(dat.xy), 
                 method.tau = "DL", hakn = FALSE)
  
  
  dir0 <- "/results/Meta/Breast"
  dir <- paste(paste(dir0, paste("Fig2_forest", paste(int[k], id[k], sep="_"), sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    25, #length
    15, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p<- forest.meta(res,
                  leftcols = c("studlab", "TE", "seTE", "lower", "upper", "pval"),
                  leftlabs = c(paste(int[k], id[k], sep="/"), "Effect", "SE", "95% CI \n Lower",
                               "95% CI \n Upper", "P value"),
                  xlab = "effect estimate", 
                  lab.e = "Intervention",
                  sortvar = TE, 
                  # pooled.totals = TRUE,
                  smlab = " ",
                  text.random = "Random effect",
                  print.I2.ci = FALSE,
                  print.Q = TRUE,
                  print.pval.Q = TRUE,
                  digits.sd = 2,
                  print.I2 = TRUE,
                  print.tau2 = TRUE,
                  comb.random =TRUE,
                  comb.fixed = TRUE,
                  text.fixed.w = TRUE,
                  colgap.forest.left = "0.5cm",
                  layout = "RevMan5", 
                  test.overall.random = TRUE,
                  test.overall.fixed = TRUE,
                  xlim = "symmetric",
                  col.square = "grey70",
                  col.inside = "grey70",
                  col.square.lines = "grey30",
                  col.diamond.random = "blue2",
                  col.diamond.fixed  = "red2",
                  ff.xlab = "bold",
                  fontsize = 11,
                  fs.heading = 11.5,
                  squaresize = 0.55,
                  scientific.pval = TRUE,
                  lty.random = NULL,
                  lty.fixed  = NULL)
  
  
  print(p)
  
  dev.off()
  
}

####################################################################
## forest plot: pan data
####################################################################
load("/results/LM/Pan/lm.pan.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Pan/lm.pan.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Pan/lm.pan.ctrpv.RData")
res.ctrp <- lm.res
load("/results/LM/Pan/lm.pan.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Pan/lm.pan.gdsc2.RData")
res.gdsc2 <- lm.res

int <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id <- c("EGFR", "ERBB2", "S100A1")

for(k in 1:length(int)){
  
  sub.res.ccle <- res.ccle[res.ccle$gene.name == id[k] & 
                             res.ccle$drug == int[k], ]
  sub.res.gcsi <- res.gcsi[res.gcsi$gene.name == id[k] & 
                             res.gcsi$drug == int[k], ]
  sub.res.gdsc1 <- res.gdsc1[res.gdsc1$gene.name == id[k] & 
                               res.gdsc1$drug == int[k], ]
  sub.res.gdsc2 <- res.gdsc2[res.gdsc2$gene.name == id[k] & 
                               res.gdsc2$drug == int[k], ]
  sub.res.ctrp <- res.ctrp[res.ctrp$gene.name == id[k] & 
                             res.ctrp$drug == int[k], ]
  
  res <- rbind(sub.res.ccle, sub.res.ctrp, sub.res.gcsi,
               sub.res.gdsc1, sub.res.gdsc2)
  study.id <- c("CCLE", "CTRP", "gCSI", "GDSC1","GDSC2")
  df.id <- paste("n", res$df + 8, sep=" = ")
  rownames(res) <- paste(study.id, df.id, sep=", ")
  
  x <-  res$estimate
  y <-  res$std
  dat.xy <- data.frame(x, y)
  rownames(dat.xy) <- rownames(res)
  res <- metagen(x, y, studlab = rownames(dat.xy), 
                 method.tau = "DL", hakn = FALSE)
  
  
  dir0 <- "/results/Meta/Pan"
  dir <- paste(paste(dir0, paste("Fig2_forest", paste(int[k], id[k], sep="_"), sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    25, #length
    15, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p<- forest.meta(res,
                  leftcols = c("studlab", "TE", "seTE", "lower", "upper", "pval"),
                  leftlabs = c(paste(int[k], id[k], sep="/"), "Effect", "SE", "95% CI \n Lower",
                               "95% CI \n Upper", "P value"),
                  xlab = "effect estimate", 
                  lab.e = "Intervention",
                  sortvar = TE, 
                  # pooled.totals = TRUE,
                  smlab = " ",
                  text.random = "Random effect",
                  print.I2.ci = FALSE,
                  print.Q = TRUE,
                  print.pval.Q = TRUE,
                  digits.sd = 2,
                  print.I2 = TRUE,
                  print.tau2 = TRUE,
                  comb.random =TRUE,
                  comb.fixed = TRUE,
                  text.fixed.w = TRUE,
                  colgap.forest.left = "0.5cm",
                  layout = "RevMan5", 
                  test.overall.random = TRUE,
                  test.overall.fixed = TRUE,
                  xlim = "symmetric",
                  col.square = "grey70",
                  col.inside = "grey70",
                  col.square.lines = "grey30",
                  col.diamond.random = "blue2",
                  col.diamond.fixed  = "red2",
                  ff.xlab = "bold",
                  fontsize = 11,
                  fs.heading = 11.5,
                  squaresize = 0.55,
                  scientific.pval = TRUE,
                  lty.random = NULL,
                  lty.fixed  = NULL)
  
  
  print(p)
  
  dev.off()
  
}


####################################################################
## volcano plot: breast 
####################################################################

load("/results/Meta/Breast/meta.breast.RData")

id.drug <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id.gene <- c("EGFR", "ERBB2", "S100A1")

for(i in 1:length(id.drug)){
  
  res <- meta.effect.size.res[meta.effect.size.res$drug == id.drug[i] &
                                meta.effect.size.res$method.meta == "DL", ]
  
  res$diffexpressed <- "NO"
  res$diffexpressed[res$TE.random > 0 & res$padj < 0.05] <- "FDR < 0.05, overall effect > 0"
  res$diffexpressed[res$TE.random < (0) & res$padj < 0.05] <- "FDR < 0.05, overall effect < 0"
  
  
  mycolors <- c( "yellow3","deepskyblue1", "grey50")
  names(mycolors) <- c("FDR < 0.05, overall effect > 0", 
                       "FDR < 0.05, overall effect < 0", 
                       "NO")
  
  res$delabel <- NA
  mod.res <- res[abs(res$TE.random) > 0 & res$padj < 0.05, ]
  mod.res <- mod.res[order(mod.res$pval.random, decreasing = FALSE),]
  id <- c(mod.res$gene.name[1:10], mod.res[mod.res$TE.random < 0 & 
                                             mod.res$pval.random < 1e-5, "gene.name"][1:10], id.gene[i])
  
  for(j in 1:length(id)){
    
    k <- which(res$gene.name == id[j])
    res$delabel[k] <- res[k, ]$gene.name
    
  }
  
  
  dir0 <- "/results/Meta/Breast"
  dir <- paste(paste(dir0, paste("Fig2_volcano", paste(id.drug[i], id.gene[i], sep="_"), sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    15, #length
    12, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <- ggplot(data=res, aes(x=TE.random, y=-log10(pval.random), 
                            col=diffexpressed)) + 
    geom_point(size = 2) + theme_minimal() +
    ylab(paste("-log10 P value", id.drug[i], sep=", ")) +  
    xlab("overall effect") +
    scale_colour_manual(values = mycolors) +
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=12,face="bold"),
          axis.text.y=element_text(size=10, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="bottom",
          legend.text = element_text(size = 8, face="bold"),
          legend.title = element_blank()) +
    geom_text_repel(aes(label=delabel),
                    size = 2,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both", 
                    seed = 2356,
                    fontface= "bold")
  
  print(p)
  
  dev.off()
  
  
}


####################################################################
## volcano plot: pan
####################################################################

load("/results/Meta/Pan/meta.pan.RData")

id.drug <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id.gene <- c("EGFR", "ERBB2", "S100A1")

for(i in 1:length(id.drug)){
  
  res <- meta.effect.size.res[meta.effect.size.res$drug == id.drug[i] &
                                meta.effect.size.res$method.meta == "DL", ]
  
  res$diffexpressed <- "NO"
  res$diffexpressed[res$TE.random > 0 & res$padj < 0.05] <- "FDR < 0.05, overall effect > 0"
  res$diffexpressed[res$TE.random < (0) & res$padj < 0.05] <- "FDR < 0.05, overall effect < 0"
  
  
  mycolors <- c( "yellow3","deepskyblue1", "grey50")
  names(mycolors) <- c("FDR < 0.05, overall effect > 0", 
                       "FDR < 0.05, overall effect < 0", 
                       "NO")
  
  res$delabel <- NA
  mod.res <- res[abs(res$TE.random) > 0 & res$padj < 0.05, ]
  mod.res <- mod.res[order(mod.res$pval.random, decreasing = FALSE),]
  id <- c(mod.res$gene.name[1:10], mod.res[mod.res$TE.random < 0 & 
                                             mod.res$pval.random < 1e-5, "gene.name"][1:10], id.gene[i])
  
  for(j in 1:length(id)){
    
    k <- which(res$gene.name == id[j])
    res$delabel[k] <- res[k, ]$gene.name
    
  }
  
  
  dir0 <- "/results/Meta/Pan"
  dir <- paste(paste(dir0, paste("Fig2_volcano", paste(id.drug[i], id.gene[i], sep="_"), sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    15, #length
    12, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <- ggplot(data=res, aes(x=TE.random, y=-log10(pval.random), 
                            col=diffexpressed)) + 
    geom_point(size = 2) + theme_minimal() +
    ylab(paste("-log10 P value", id.drug[i], sep=", ")) +  
    xlab("overall effect") +
    scale_colour_manual(values = mycolors) +
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=12,face="bold"),
          axis.text.y=element_text(size=10, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="bottom",
          legend.text = element_text(size = 8, face="bold"),
          legend.title = element_blank()) +
    geom_text_repel(aes(label=delabel),
                    size = 2,
                    color = "black",
                    min.segment.length = 0,
                    na.rm = TRUE, direction = "both", 
                    seed = 2356,
                    fontface= "bold")
  
  print(p)
  
  dev.off()
  
  
}

################################################################
## Jaccard Index: breast 
################################################################

load("/results/Meta/Breast/meta.breast.RData")

k<- 100 
id.drug <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(j in 1:length(id.drug)){
  
  res.fisher <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.fisher <- res.fisher[order(res.fisher$fisher.p), ]
  sub.fisher <- res.fisher[1:k, c("gene.name", "drug")]
  sub.fisher <- data.frame(method = "Fisher", sub.fisher)
  
  res.stouffer <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.stouffer <- res.stouffer[order(res.stouffer$stouffer.p), ]
  sub.stouffer <- res.stouffer[1:k, c("gene.name", "drug")]
  sub.stouffer <- data.frame(method = "Stouffer", sub.stouffer)
  
  res.SJ <- meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                   meta.effect.size.res$method.meta == "SJ", ]
  res.SJ <- res.SJ[order(res.SJ$pval.random), ]
  sub.SJ <- res.SJ[1:k, c("gene.name", "drug")]
  sub.SJ <- data.frame(method = "SJ", sub.SJ)
  
  res.DL <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "DL", ]
  res.DL <- res.DL[order(res.DL$pval.random), ]
  sub.DL <- res.DL[1:k, c("gene.name", "drug")]
  sub.DL <- data.frame(method = "DL", sub.DL)
  
  res.PM <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "PM", ]
  res.PM <- res.PM[order(res.PM$pval.random), ]
  sub.PM <- res.PM[1:k, c("gene.name", "drug")]
  sub.PM <- data.frame(method = "PM", sub.PM)
  
  res.HS <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "HS", ]
  res.HS <- res.HS[order(res.HS$pval.random), ]
  sub.HS <- res.HS[1:k, c("gene.name", "drug")]
  sub.HS <- data.frame(method = "HS", sub.HS)
  
  res.HE <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "HE", ]
  res.HE <- res.HE[order(res.HE$pval.random), ]
  sub.HE <- res.HE[1:k, c("gene.name", "drug")]
  sub.HE <- data.frame(method = "HE", sub.HE)
  
  res.EB <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "EB", ]
  res.EB <- res.EB[order(res.EB$pval.random), ]
  sub.EB <- res.EB[1:k, c("gene.name", "drug")]
  sub.EB <- data.frame(method = "EB", sub.EB)
  
  data <- rbind(sub.fisher, sub.stouffer, sub.SJ, 
                sub.DL, sub.PM, sub.HS, sub.HE, sub.EB) 
  
  id <- unique(data$method)
  
  sub.data <-  lapply(1:length(id), function(i){
    
    
    dat <- data[data$method == id[i], ]
    
    mod.dat <- sapply(1:length(id), function(k){
      
      jaccard(data[data$method == id[k], "gene.name"],
              dat[, "gene.name"])
      
      
    })
    
    mod.dat
    
  })
  
  jc.index <- do.call(rbind, sub.data)
  colnames(jc.index) <- id
  rownames(jc.index) <- id
  
  dir0 <- "/results/Meta/Breast"
  dir <- paste(paste(dir0, paste("SFig4_JaccardIndex_cont", id.drug[j], sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    16, #length
    17, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  col <- colorRampPalette(c("deepskyblue1", "white","deepskyblue4"))
  par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
  plot(as.assoc(jc.index), col=col, 
       breaks= range(jc.index), 
       cex=0.60, key = NULL, 
       main = id.drug[j], xlab = "", ylab = "",
       cex.axis = 0.68)
  
  dev.off()
  
}


####################################################################
## Jaccard Index: Pan data 
####################################################################

load("/results/Meta/Pan/meta.pan.RData")

k<- 100 
id.drug <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(j in 1:length(id.drug)){
  
  res.fisher <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.fisher <- res.fisher[order(res.fisher$fisher.p), ]
  sub.fisher <- res.fisher[1:k, c("gene.name", "drug")]
  sub.fisher <- data.frame(method = "Fisher", sub.fisher)
  
  res.stouffer <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.stouffer <- res.stouffer[order(res.stouffer$stouffer.p), ]
  sub.stouffer <- res.stouffer[1:k, c("gene.name", "drug")]
  sub.stouffer <- data.frame(method = "Stouffer", sub.stouffer)
  
  res.SJ <- meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                   meta.effect.size.res$method.meta == "SJ", ]
  res.SJ <- res.SJ[order(res.SJ$pval.random), ]
  sub.SJ <- res.SJ[1:k, c("gene.name", "drug")]
  sub.SJ <- data.frame(method = "SJ", sub.SJ)
  
  res.DL <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "DL", ]
  res.DL <- res.DL[order(res.DL$pval.random), ]
  sub.DL <- res.DL[1:k, c("gene.name", "drug")]
  sub.DL <- data.frame(method = "DL", sub.DL)
  
  res.PM <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "PM", ]
  res.PM <- res.PM[order(res.PM$pval.random), ]
  sub.PM <- res.PM[1:k, c("gene.name", "drug")]
  sub.PM <- data.frame(method = "PM", sub.PM)
  
  res.HS <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "HS", ]
  res.HS <- res.HS[order(res.HS$pval.random), ]
  sub.HS <- res.HS[1:k, c("gene.name", "drug")]
  sub.HS <- data.frame(method = "HS", sub.HS)
  
  res.HE <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "HE", ]
  res.HE <- res.HE[order(res.HE$pval.random), ]
  sub.HE <- res.HE[1:k, c("gene.name", "drug")]
  sub.HE <- data.frame(method = "HE", sub.HE)
  
  res.EB <-  meta.effect.size.res[meta.effect.size.res$drug == id.drug[j] &
                                    meta.effect.size.res$method.meta == "EB", ]
  res.EB <- res.EB[order(res.EB$pval.random), ]
  sub.EB <- res.EB[1:k, c("gene.name", "drug")]
  sub.EB <- data.frame(method = "EB", sub.EB)
  
  data <- rbind(sub.fisher, sub.stouffer, sub.SJ, 
                sub.DL, sub.PM, sub.HS, sub.HE, sub.EB) 
  
  id <- unique(data$method)
  
  sub.data <-  lapply(1:length(id), function(i){
    
    
    dat <- data[data$method == id[i], ]
    
    mod.dat <- sapply(1:length(id), function(k){
      
      jaccard(data[data$method == id[k], "gene.name"],
              dat[, "gene.name"])
      
      
    })
    
    mod.dat
    
  })
  
  jc.index <- do.call(rbind, sub.data)
  colnames(jc.index) <- id
  rownames(jc.index) <- id
  
  dir0 <- "/results/Meta/Pan"
  dir <- paste(paste(dir0, paste("SFig4_JaccardIndex_cont", id.drug[j], sep="_"), sep="/"), "jpeg", sep=".")
  
  Cairo::Cairo(
    16, #length
    17, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  col <- colorRampPalette(c("deepskyblue1", "white","deepskyblue4"))
  par(mar=c(5.1, 4.1, 4.1, 4.1)) # adapt margins
  plot(as.assoc(jc.index), col=col, 
       breaks= range(jc.index), 
       cex=0.60, key = NULL, 
       main = id.drug[j], xlab = "", ylab = "",
       cex.axis = 0.68)
  
  dev.off()
  
}

####################################################################
## Breast data: I2 and DE results  
####################################################################

load("/results/Meta/Breast/meta.breast.RData")

id.drug <- c("Erlotinib", "Lapatinib", "Paclitaxel")

res.perc <- lapply(1:length(id.drug), function(k){
  
  res <- meta.effect.size.res[meta.effect.size.res$drug == id.drug[k] &
                                meta.effect.size.res$method.meta == "DL"  , ]
  res$I2.id <- ifelse(res$I2 <= 0.5 & res$pval.Q > 0.1, 
                      "not substantial", "substantial")
  res$sig.id <- ifelse(res$padj < 0.05, "DE", "non-DE")
  perc.val <- round(percentages(table(res$I2.id, res$sig.id)))
  
  dat <- data.frame(sig = c("DE", "DE", "non-DE", "non-DE"),
                    I2 = c("not substantial", "substantial", "not substantial", "substantial"),
                    y = c(perc.val[1,1], perc.val[2,1], perc.val[1,2], perc.val[2,2]),
                    drug = id.drug[k])
  dat
  
})


dat <- do.call(rbind, res.perc)
dat$I2 <- factor(dat$I2)
levels(dat$I2) <- c("I2 < 50%", "I2 > 50%")
dat$sig <- factor(dat$sig)
levels(dat$sig) <- c("FDR < 0.05", "FDR > 0.05")

dir0 <- "/results/Meta/Breast"
dir <- paste(paste(dir0, "SFig6_I2_DE", sep="/"), "jpg", sep=".")

Cairo::Cairo(
  15, #length
  10, #width
  file = dir,
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

p <-ggplot(dat , aes(y =  y, x = I2, fill = sig)) +
  geom_bar(width = 0.3, stat="identity", colour = "grey") +
  facet_wrap(drug ~ .) +
  scale_fill_manual(values=c("turquoise3", "ivory4")) +
  xlab("") +
  ylab("percentage of genes") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_blank())

print(p)

dev.off()


####################################################################
## Pan data: I2 and DE results  
####################################################################

load("/results/Meta/Pan/meta.pan.RData")

id.drug <- c("Erlotinib", "Lapatinib", "Paclitaxel")

res.perc <- lapply(1:length(id.drug), function(k){
  
  res <- meta.effect.size.res[meta.effect.size.res$drug == id.drug[k] &
                                meta.effect.size.res$method.meta == "DL"  , ]
  res$I2.id <- ifelse(res$I2 <= 0.5 & res$pval.Q > 0.1, 
                      "not substantial", "substantial")
  res$sig.id <- ifelse(res$padj < 0.05, "DE", "non-DE")
  perc.val <- round(percentages(table(res$I2.id, res$sig.id)))
  
  dat <- data.frame(sig = c("DE", "DE", "non-DE", "non-DE"),
                    I2 = c("not substantial", "substantial", "not substantial", "substantial"),
                    y = c(perc.val[1,1], perc.val[2,1], perc.val[1,2], perc.val[2,2]),
                    drug = id.drug[k])
  dat
  
})


dat <- do.call(rbind, res.perc)
dat$I2 <- factor(dat$I2)
levels(dat$I2) <- c("I2 < 50%", "I2 > 50%")
dat$sig <- factor(dat$sig)
levels(dat$sig) <- c("FDR < 0.05", "FDR > 0.05")

dir0 <- "/results/Meta/Pan"
dir <- paste(paste(dir0, "SFig6_I2_DE", sep="/"), "jpg", sep=".")

Cairo::Cairo(
  15, #length
  10, #width
  file = dir,
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

p <-ggplot(dat , aes(y =  y, x = I2, fill = sig)) +
  geom_bar(width = 0.3, stat="identity", colour = "grey") +
  facet_wrap(drug ~ .) +
  scale_fill_manual(values=c("turquoise3", "ivory4")) +
  #ggtitle(plot.id) + 
  xlab("") +
  ylab("percentage of genes") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=10,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 10, face="bold"),
        legend.title = element_blank())

print(p)

dev.off()

########################################################################################
########################################################################################
################################### meta-analysis: Bayesian analyses ###################
########################################################################################
########################################################################################

#########################################################
##  forest plot breast: top genes and common drugs
#########################################################

load("/results/Meta/Breast/meta.breast.RData")
meta.dl <- meta.effect.size.res
load("/results/Meta/Breast/meta.breast.bayes.RData")
meta.bayes <- meta.effect.size.res

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

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
gene.id <- c("EGFR", "ERBB2", "S100A1")

lm.dat <- list(res.ccle, res.ctrpv, res.gcsi, res.gray, res.gdsc1, res.gdsc2,  res.uhn)
study <- c("CCLE", "CTRP", "gCSI", "GRAY", "GDSC1", "GDSC2", "UHN")

combine.lm.res <- lapply(1:length(drugs), function(k){
  
  res <- lapply(1:length(lm.dat), function(i){
    
    data.frame(study.id = study[i], lm.dat[[i]][lm.dat[[i]]$drug == drugs[k], ]) 
    
  })
  
  do.call(rbind, res)
  
})

combine.lm.res <- do.call(rbind, combine.lm.res)

for(i in 1:length(drugs)){
  
  sub.res.lm <- combine.lm.res[combine.lm.res$gene.name == gene.id[i] & combine.lm.res$drug == drugs[i], ]
  res.dl <- meta.dl[meta.dl$drug == drugs[i] & 
                      meta.dl$gene.name == gene.id[i] &
                      meta.dl$method.meta == "DL", ]
  
  res.bayes <- meta.bayes[meta.bayes$drug == drugs[i] & meta.bayes$gene.name == gene.id[i], ]
  
  coef.res <- c(sub.res.lm$estimate,
                res.dl$TE.random, 
                res.bayes$coef)
  
  
  ci.low <- c(sub.res.lm$ci.low,
              res.dl$lower.random,
              res.bayes$ci.low)
  
  ci.up <- c(sub.res.lm$ci.high,
             res.dl$upper.random, 
             res.bayes$ci.up)
  
  ci.lower <- paste("[", round(ci.low,2), sep="")
  ci.upper <- paste(round(ci.up,2), "]",sep="")
  ci <- paste(ci.lower, ci.upper, sep=",")
  beta.ci <- paste(round(coef.res, 2), 
                   ci, sep=" ")
  
  
  cochrane_from_rmeta <- 
    structure(list(
      mean  = c(NA, round(coef.res,2)), 
      lower = c(NA, round(ci.low,2)),
      upper = c(NA, round(ci.up, 2))),
      .Names = c("Beta", "lower", "upper"), 
      row.names = c(NA, -10L), 
      class = "data.frame")
  
  i2.dl <- paste(round(res.dl$I2, 2) * 100, "%", sep = "")
  i2.dl <- paste("I2", i2.dl, sep="=")
  i2.bayes <- paste(res.bayes$I2, "%", sep = "")
  i2.bayes <- paste("I2", i2.bayes, sep="=")
  
  tabletext<-cbind(
    c(paste(drugs[i], gene.id[i], sep="/"), c(paste(study, "n=19", sep = ", "),
                                              paste("DerSimonian-Laird (DL)", i2.dl, sep=", "),  
                                              paste("Jeffreys", i2.bayes, sep=", "))),
    c("Beta [95% CI]", beta.ci))
  
  
  dir0 <- "/results/MetaVisualization/Breast"
  dir <- paste(dir0, paste(paste("SFig7", drugs[i], sep = "_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    20, #length
    15, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <- forestplot(tabletext, 
                  hrzl_lines = gpar(col="grey50", lwd=2),
                  cochrane_from_rmeta, new_page = TRUE,
                  is.summary=c(TRUE,rep(FALSE,7), TRUE, TRUE),
                  #clip=c(0.1,1.1), 
                  xlog=FALSE,
                  xticks.digits = 1,
                  boxsize = .3, # We set the box size to better visualize the type
                  line.margin = .2, # We need to add this to avoid crowding
                  xlab = "effect estimate",
                  txt_gp = fpTxtGp(xlab  = gpar(cex = 0.9, fontface = "bold"),
                                   label  =  list(gpar(cex = 0.9, fontface="bold", col = "grey20"),
                                                  gpar(cex = 0.8),
                                                  gpar(cex = 0.8))),
                  col=fpColors(box="grey40", line="grey10",
                               summary ="blue" ))
  
  
  print(p)
  
  
  dev.off()
  
  
}

#########################################################
##  forest plot pan: top genes and common drugs
#########################################################

load("/results/Meta/Pan/meta.pan.RData")
meta.dl <- meta.effect.size.res
load("/results/Meta/Pan/meta.pan.bayes.RData")
meta.bayes <- meta.effect.size.res

load("/results/LM/Pan/lm.pan.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Pan/lm.pan.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Pan/lm.pan.ctrpv.RData")
res.ctrpv <- lm.res
load("/results/LM/Pan/lm.pan.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Pan/lm.pan.gdsc2.RData")
res.gdsc2 <- lm.res

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
gene.id <- c("EGFR", "ERBB2", "S100A1")

lm.dat <- list(res.ccle, res.ctrpv, res.gcsi,  res.gdsc1, res.gdsc2)
study <- c("CCLE", "CTRP", "gCSI", "GDSC1", "GDSC2")

combine.lm.res <- lapply(1:length(drugs), function(k){
  
  res <- lapply(1:length(lm.dat), function(i){
    
    data.frame(study.id = study[i], lm.dat[[i]][lm.dat[[i]]$drug == drugs[k], ]) 
    
  })
  
  do.call(rbind, res)
  
})

combine.lm.res <- do.call(rbind, combine.lm.res)

for(i in 1:length(drugs)){
  
  sub.res.lm <- combine.lm.res[combine.lm.res$gene.name == gene.id[i] & combine.lm.res$drug == drugs[i], ]
  res.dl <- meta.dl[meta.dl$drug == drugs[i] & 
                      meta.dl$gene.name == gene.id[i] &
                      meta.dl$method.meta == "DL", ]
  
  res.bayes <- meta.bayes[meta.bayes$drug == drugs[i] & meta.bayes$gene.name == gene.id[i], ]
  
  coef.res <- c(sub.res.lm$estimate,
                res.dl$TE.random, 
                res.bayes$coef)
  
  ci.low <- c(sub.res.lm$ci.low,
              res.dl$lower.random,
              res.bayes$ci.low)
  
  ci.up <- c(sub.res.lm$ci.high,
             res.dl$upper.random, 
             res.bayes$ci.up)
  
  ci.lower <- paste("[", round(ci.low,2), sep="")
  ci.upper <- paste(round(ci.up,2), "]",sep="")
  ci <- paste(ci.lower, ci.upper, sep=",")
  beta.ci <- paste(round(coef.res, 2), 
                   ci, sep=" ")
  
  cochrane_from_rmeta <- 
    structure(list(
      mean  = c(NA, round(coef.res,2)), 
      lower = c(NA, round(ci.low,2)),
      upper = c(NA, round(ci.up, 2))),
      .Names = c("Beta", "lower", "upper"), 
      row.names = c(NA, -8L), 
      class = "data.frame")
  
  i2.dl <- paste(round(res.dl$I2, 2) * 100, "%", sep = "")
  i2.dl <- paste("I2", i2.dl, sep="=")
  i2.bayes <- paste(res.bayes$I2, "%", sep = "")
  i2.bayes <- paste("I2", i2.bayes, sep="=")
  
  tabletext<-cbind(
    c(paste(drugs[i], gene.id[i], sep="/"), c(paste(study, "n=168", sep = ", "),
                                              paste("DerSimonian-Laird (DL)", i2.dl, sep=", "),  
                                              paste("Jeffreys", i2.bayes, sep=", "))),
    c("Beta [95% CI]", beta.ci))
  
  dir0 <- "/results/MetaVisualization/Pan"
  dir <- paste(dir0, paste(paste("SFig7", drugs[i], sep = "_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    20, #length
    15, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <- forestplot(tabletext, 
                  hrzl_lines = gpar(col="grey50", lwd=2),
                  cochrane_from_rmeta, new_page = TRUE,
                  is.summary=c(TRUE,rep(FALSE,5), TRUE, TRUE),
                  #clip=c(0.1,1.1), 
                  xlog=FALSE,
                  xticks.digits = 1,
                  boxsize = .3, # We set the box size to better visualize the type
                  line.margin = .2, # We need to add this to avoid crowding
                  xlab = "effect estimate",
                  txt_gp = fpTxtGp(xlab  = gpar(cex = 0.9, fontface = "bold"),
                                   label  =  list(gpar(cex = 0.9, fontface="bold", col = "grey20"),
                                                  gpar(cex = 0.8),
                                                  gpar(cex = 0.8))),
                  col=fpColors(box="grey40", line="grey10",
                               summary ="blue" ))
  
  
  print(p)
  
  
  dev.off()
  
  
}

########################################################################################
## Breast data: compare Bayesian and DL meta-analyses (I2 and overall effect size) 
########################################################################################

load("/results/Meta/Breast/meta.breast.RData")
meta.dl <- meta.effect.size.res
load("/results/Meta/Breast/meta.breast.bayes.RData")
meta.bayes <- meta.effect.size.res

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(i in 1:length(drugs)){
  
  res.bayes <- meta.bayes[meta.bayes$drug == drugs[i], ]
  res.bayes <- res.bayes[order(res.bayes$gene.name), ]
  
  res.dl <- meta.dl[meta.dl$drug == drugs[i] &
                      meta.dl$method.meta == "DL", ]
  res.dl <- res.dl[order(res.dl$gene.name), ]
  
  ci.dl <- res.dl$upper.random - res.dl$lower.random
  ci.bayes <- res.bayes$ci.up - res.bayes$ci.low
  
  i2.dl <- res.dl$I2
  i2.bayes <- res.bayes$I2/100
  
  dat <- data.frame(I2 = c(i2.dl, i2.bayes),
                    ci = c(ci.dl, ci.bayes),
                    group = c(rep("DerSimonian-Laird (DL)", length(ci.dl)),
                              rep("Jeffreys", length(ci.bayes))))
  
  
  dir0 <- "/results/MetaVisualization/Breast"
  dir <- paste(dir0, paste(paste("SFig8_ci", drugs[i], sep = "_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    11, #length
    11, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <-ggplot(dat , aes(group, ci, fill = group)) +
    geom_violin(width = 1) +
    stat_summary(fun.y=median, geom="point", size=2, color="black") +
    scale_fill_manual(values=c("maroon", "seashell3")) +
    #geom_point(position = position_jitter(width = 0.2, seed=353), size=1.4) +
    #ggtitle("KW P value = 0.998") + 
    xlab("") +
    ylab(paste("length of 95% credible or confidence interval", drugs[i], sep=", ") )+
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=8,face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
  dir0 <- "/results/MetaVisualization/Breast"
  dir <- paste(dir0, paste(paste("SFig8_i2", drugs[i], sep = "_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    11, #length
    11, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <-ggplot(dat , aes(group, I2, fill = group)) +
    geom_violin(width = 1) +
    stat_summary(fun.y=median, geom="point", size=2, color="black") +
    scale_fill_manual(values=c("maroon", "seashell3")) +
    xlab("") +
    ylab(paste("heterogeneity estimate", drugs[i], sep=", ")) +
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=8,face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
}


####################################################################################
## Pan data: compare Bayesian and DL meta-analyses (I2 and overall effect size) 
####################################################################################

load("/results/Meta/Pan/meta.pan.RData")
meta.dl <- meta.effect.size.res
load("/results/Meta/Pan/meta.pan.bayes.RData")
meta.bayes <- meta.effect.size.res

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(i in 1:length(drugs)){
  
  res.bayes <- meta.bayes[meta.bayes$drug == drugs[i], ]
  res.bayes <- res.bayes[order(res.bayes$gene.name), ]
  
  res.dl <- meta.dl[meta.dl$drug == drugs[i] &
                      meta.dl$method.meta == "DL", ]
  res.dl <- res.dl[order(res.dl$gene.name), ]
  
  
  ci.dl <- res.dl$upper.random - res.dl$lower.random
  ci.bayes <- res.bayes$ci.up - res.bayes$ci.low
  
  i2.dl <- res.dl$I2
  i2.bayes <- res.bayes$I2/100
  
  dat <- data.frame(I2 = c(i2.dl, i2.bayes),
                    ci = c(ci.dl, ci.bayes),
                    group = c(rep("DerSimonian-Laird (DL)", length(ci.dl)),
                              rep("Jeffreys", length(ci.bayes))))
  
  
  dir0 <- "/results/MetaVisualization/Pan"
  dir <- paste(dir0, paste(paste("SFig8_ci", drugs[i], sep = "_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    11, #length
    11, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <-ggplot(dat , aes(group, ci, fill = group)) +
    geom_violin(width = 1) +
    stat_summary(fun.y=median, geom="point", size=2, color="black") +
    scale_fill_manual(values=c("maroon", "seashell3")) +
    #geom_point(position = position_jitter(width = 0.2, seed=353), size=1.4) +
    #ggtitle("KW P value = 0.998") + 
    xlab("") +
    ylab(paste("length of 95% credible or confidence interval", drugs[i], sep=", ") )+
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=8,face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
  dir0 <- "/results/MetaVisualization/Pan"
  dir <- paste(dir0, paste(paste("SFig8_i2", drugs[i], sep = "_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    11, #length
    11, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
  
  p <-ggplot(dat , aes(group, I2, fill = group)) +
    geom_violin(width = 1) +
    stat_summary(fun.y=median, geom="point", size=2, color="black") +
    scale_fill_manual(values=c("maroon", "seashell3")) +
    xlab("") +
    ylab(paste("heterogeneity estimate", drugs[i], sep=", ")) +
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=8,face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
}


