source("/code/Prog_LinearModel_Meta_fun.R")
source("/code/Prog_correlation.run.R")
source("/code/Prog_upset_meta.R")

dir.create("/results/LM")
dir.create("/results/LM/Breast")

dir.create("/results/MetaNoIndept")
dir.create("/results/MetaNoIndept/Breast")

################################################################################################
################################################################################################
#################################### linear model fitting ######################################
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

#save(meta.effect.size.res, meta.pvalue.res, 
#     file = "/results/Meta/Breast/meta.breast.RData")

################################################################################################
################################################################################################
#################################### non-independent analyses ##################################
################################################################################################
################################################################################################

load("/data/breast/DrugMeta.breast.zfixsample.RData")
load("/data/breast/DrugMeta.breast.fixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

# NOTE: change the loop to length(id.gene)
id.gene <- z.exprs.gcsi.data$gene_name

set.seed(135)
id.gene <- c(sample(id.gene, 50),
             "CHMP7", "FBXO3", "DUOX2", "EGFR", "CD63", "ZNF143",
             "PSMB3", "C3orf52", "PITX3", "ERBB2", "CD63", "CD84",
             "IGF2BP3", "TGM3", "ABCG4", "S100A1", "COL11A1", "SCRT1")


dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                            meta.effect.size.res$drug == drugs[l] &
                                            meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               tau =  meta.res.indept$tau2,
                               I2 = meta.res.indept$I2,
                               pval.Q = meta.res.indept$pval.Q,
                               te.val.r = meta.res.indept$TE.random,
                               sete.val.r = meta.res.indept$seTE.random,
                               te.val.f = meta.res.indept$TE.fixed,
                               sete.val.f = meta.res.indept$seTE.fixed,
                               pval.r = meta.res.indept$pval.random,
                               pval.f = meta.res.indept$pval.fixed)
    
    ##############################
    ## One duplicate analysis
    ##############################
    
    res.all <- lapply(1:length(study.id), function(m){
      
      id.missing <- study.id[m]
      id.known <- study.id[study.id != id.missing]
      mod.dat.exprs <- dat.exprs
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        mod.dat.exprs
        
      })
      
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        res <- data.frame(study.unknown = id.missing, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- metagen(x, y, studlab = rownames(dat.xy), 
                       method.tau = "DL")
        
        res
        
      })
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   tau =  meta.res[[k]]$tau2,
                   I2 = meta.res[[k]]$I2,
                   pval.Q = meta.res[[k]]$pval.Q,
                   te.val.r = meta.res[[k]]$TE.random,
                   sete.val.r = meta.res[[k]]$seTE.random,
                   te.val.f = meta.res[[k]]$TE.fixed,
                   sete.val.f = meta.res[[k]]$seTE.fixed,
                   pval.r = meta.res[[k]]$pval.random,
                   pval.f = meta.res[[k]]$pval.fixed)
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/one.dup.RData")

################################################################
################################################################
## two duplication
################################################################
#################################################################


dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    # load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               tau =  meta.res.indept$tau2,
                               I2 = meta.res.indept$I2,
                               pval.Q = meta.res.indept$pval.Q,
                               te.val.r = meta.res.indept$TE.random,
                               sete.val.r = meta.res.indept$seTE.random,
                               te.val.f = meta.res.indept$TE.fixed,
                               sete.val.f = meta.res.indept$seTE.fixed,
                               pval.r = meta.res.indept$pval.random,
                               pval.f = meta.res.indept$pval.fixed)
    
  ####################################################################
  ## two duplicate analysis
####################################################################
  
  id <- lapply(1:(length(study.id)-1), function(k){
    
    id1 <- study.id[k]
    
    
    id.dup <-  sapply((k+1):length(study.id), function(i){
      
      cbind(id1, study.id[i])
      
    })
    
    id.dup
    
  })
  
  
  res.all <- lapply(1:length(id), function(m){
    
    id.dup <- id[[m]]
    
    sub.res <-  lapply(1:ncol(id.dup), function(s){
      
      id.missing <- id.dup[, s]
      id.known <- study.id[!(study.id %in% id.missing)]
      mod.dat.exprs <- dat.exprs
      
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        
        mod.dat.exprs
        
      })
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        id.miss <- paste(id.missing[1], id.missing[2], sep="-")
        res <- data.frame(study.unknown = id.miss, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- metagen(x, y, studlab = rownames(dat.xy), 
                       method.tau = "DL")
        
        res
        
      }) 
      
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   tau =  meta.res[[k]]$tau2,
                   I2 = meta.res[[k]]$I2,
                   pval.Q = meta.res[[k]]$pval.Q,
                   te.val.r = meta.res[[k]]$TE.random,
                   sete.val.r = meta.res[[k]]$seTE.random,
                   te.val.f = meta.res[[k]]$TE.fixed,
                   sete.val.f = meta.res[[k]]$seTE.fixed,
                   pval.r = meta.res[[k]]$pval.random,
                   pval.f = meta.res[[k]]$pval.fixed)
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    
    
    do.call(rbind, sub.res)
    
    
  })
  
  meta.non.indept <- do.call(rbind, res.all)
  dat.gene <- rbind( meta.indept,  meta.non.indept)
  data.frame(gene_name = id.gene[t],
             drug = drugs[l], 
             dat.gene)
  
})  

  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/two.dup.RData")

#############################################################
#############################################################
## three duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
##############################
## Independent analysis
#############################
    
  #  load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               tau =  meta.res.indept$tau2,
                               I2 = meta.res.indept$I2,
                               pval.Q = meta.res.indept$pval.Q,
                               te.val.r = meta.res.indept$TE.random,
                               sete.val.r = meta.res.indept$seTE.random,
                               te.val.f = meta.res.indept$TE.fixed,
                               sete.val.f = meta.res.indept$seTE.fixed,
                               pval.r = meta.res.indept$pval.random,
                               pval.f = meta.res.indept$pval.fixed)

###################################################################
## Three duplicate analysis
####################################################################

id.res <- lapply(1:5, function(j){
  
  id.mod <- lapply(1:(length(study.id)-(j+1)), function(k){
    
    id1 <- study.id[k]
    id2 <- study.id[k+j]
    
    id.dup <-  sapply((k+(j+1)):length(study.id), function(i){
      
      cbind(id1, id2, study.id[i])
      
    })
    
    id.dup
    
  })
  
  id.mod
  
})

mod.id <- lapply(1:length(id.res), function(k){
  
  do.call(cbind, id.res[[k]])
  
})

id <- do.call(cbind, mod.id)


res.all <- lapply(1:ncol(id), function(m){
  
    id.dup <- id[,m]
    id.missing <- id.dup
    id.known <- study.id[!(study.id %in% id.missing)]
    mod.dat.exprs <- dat.exprs
   
    
    dat.exprs.dup <- lapply(1:length(id.known), function(k){
      
      mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
      
      mod.dat.exprs
      
    })
    
    lm.res <- lapply(1:length(dat.exprs.dup), function(i){
      
      dat.expr <- dat.exprs.dup[[i]]
      
      res <- lapply(1:length(study.id), function(k){
        
        y <- dat.aac[, study.id[k]]
        y <- scale(as.numeric(y), center = T, scale = T)[,1]
        x <- dat.expr[ ,study.id[k]]
        
        fit <- lm(y ~ as.numeric(x))
        coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
        r.squared <- summary(fit)$r.squared 
        df <- fit$df.residual
        coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
        rownames(coef.fit) <- NULL
        colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
        coef.fit
        
      })
      
      res <- do.call(rbind, res)
      df.id <- paste("n", res$df + 2, sep=" = ")
      rownames(res) <- paste(study.id, df.id, sep=", ")
      id.miss <- paste(id.missing[1], id.missing[2], id.missing[3], sep="-")
      res <- data.frame(study.unknown = id.miss, 
                        study.dup = id.known[i], res)
      res
    
  })
  
 
    meta.res <- lapply(1:length(lm.res), function(j){
      
      x <-  lm.res[[j]]$estimate
      y <-  lm.res[[j]]$std
      dat.xy <- data.frame(x, y)
      
      rownames(dat.xy) <- rownames(lm.res)
      
      res <- metagen(x, y, studlab = rownames(dat.xy), 
                     method.tau = "DL")
      
      res
      
    }) 
 
    
    id.missing.dup <- lapply(1:length(lm.res), function(j){
      
      data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                 id.duplicate = unique(lm.res[[j]]$study.dup))
      
    })
    
    id.miss.dup <- do.call(rbind, id.missing.dup)
    
    meta.result <- lapply(1:length(meta.res), function(k){
      
      data.frame(id.miss.dup[k,],
                 tau =  meta.res[[k]]$tau2,
                 I2 = meta.res[[k]]$I2,
                 pval.Q = meta.res[[k]]$pval.Q,
                 te.val.r = meta.res[[k]]$TE.random,
                 sete.val.r = meta.res[[k]]$seTE.random,
                 te.val.f = meta.res[[k]]$TE.fixed,
                 sete.val.f = meta.res[[k]]$seTE.fixed,
                 pval.r = meta.res[[k]]$pval.random,
                 pval.f = meta.res[[k]]$pval.fixed)
      
    })
    
    do.call(rbind, meta.result)
   
  })
  

meta.non.indept <- do.call(rbind, res.all)
dat.gene <- rbind( meta.indept,  meta.non.indept)
data.frame(gene_name = id.gene[t],
           drug = drugs[l], 
           dat.gene)

  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/three.dup.RData")

#############################################################
#############################################################
## four duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
  #  load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               tau =  meta.res.indept$tau2,
                               I2 = meta.res.indept$I2,
                               pval.Q = meta.res.indept$pval.Q,
                               te.val.r = meta.res.indept$TE.random,
                               sete.val.r = meta.res.indept$seTE.random,
                               te.val.f = meta.res.indept$TE.fixed,
                               sete.val.f = meta.res.indept$seTE.fixed,
                               pval.r = meta.res.indept$pval.random,
                               pval.f = meta.res.indept$pval.fixed)

####################################################################
## Four duplicate analysis
####################################################################

id.res <- lapply(1:4, function(j){
  
  id.mod <- lapply(1:(length(study.id)-(j+2)), function(k){
    
    id1 <- study.id[k]
    id2 <- study.id[k+j]
    id3 <- study.id[k+j+1]
    
    id.dup <-  sapply((k+(j+2)):length(study.id), function(i){
      
      cbind(id1, id2, id3, study.id[i])
      
    })
    
    id.dup
    
  })
  
  id.mod
  
})

mod.id <- lapply(1:length(id.res), function(k){
  
  do.call(cbind, id.res[[k]])
  
})

id <- do.call(cbind, mod.id)

res.all <- lapply(1:ncol(id), function(m){
  
    id.dup <- id[,m]
    id.missing <- id.dup
    id.known <- study.id[!(study.id %in% id.missing)]
    mod.dat.exprs <- dat.exprs
   
    
    dat.exprs.dup <- lapply(1:length(id.known), function(k){
      
      mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
      
      mod.dat.exprs
      
    })
    
    lm.res <- lapply(1:length(dat.exprs.dup), function(i){
      
      dat.expr <- dat.exprs.dup[[i]]
      
      res <- lapply(1:length(study.id), function(k){
        
        y <- dat.aac[, study.id[k]]
        y <- scale(as.numeric(y), center = T, scale = T)[,1]
        x <- dat.expr[ ,study.id[k]]
        
        fit <- lm(y ~ as.numeric(x))
        coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
        r.squared <- summary(fit)$r.squared 
        df <- fit$df.residual
        coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
        rownames(coef.fit) <- NULL
        colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
        coef.fit
        
      })
      
      res <- do.call(rbind, res)
      df.id <- paste("n", res$df + 2, sep=" = ")
      rownames(res) <- paste(study.id, df.id, sep=", ")
      id.miss <- paste(id.missing[1], id.missing[2], 
                       id.missing[3], id.missing[4], sep="-")
      res <- data.frame(study.unknown = id.miss, 
                        study.dup = id.known[i], res)
      res
    
  })
  
 
    meta.res <- lapply(1:length(lm.res), function(j){
      
      x <-  lm.res[[j]]$estimate
      y <-  lm.res[[j]]$std
      dat.xy <- data.frame(x, y)
      
      rownames(dat.xy) <- rownames(lm.res)
      
      res <- metagen(x, y, studlab = rownames(dat.xy), 
                     method.tau = "DL")
      
      res
      
    }) 
 
    
    id.missing.dup <- lapply(1:length(lm.res), function(j){
      
      data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                 id.duplicate = unique(lm.res[[j]]$study.dup))
      
    })
    
    id.miss.dup <- do.call(rbind, id.missing.dup)
    
    meta.result <- lapply(1:length(meta.res), function(k){
      
      data.frame(id.miss.dup[k,],
                 tau =  meta.res[[k]]$tau2,
                 I2 = meta.res[[k]]$I2,
                 pval.Q = meta.res[[k]]$pval.Q,
                 te.val.r = meta.res[[k]]$TE.random,
                 sete.val.r = meta.res[[k]]$seTE.random,
                 te.val.f = meta.res[[k]]$TE.fixed,
                 sete.val.f = meta.res[[k]]$seTE.fixed,
                 pval.r = meta.res[[k]]$pval.random,
                 pval.f = meta.res[[k]]$pval.fixed)
      
    })
    
    do.call(rbind, meta.result)
   
  })
  


meta.non.indept <- do.call(rbind, res.all)
dat.gene <- rbind( meta.indept,  meta.non.indept)
data.frame(gene_name = id.gene[t],
           drug = drugs[l], 
           dat.gene)

  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/four.dup.RData")

#############################################################
#############################################################
## five duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
  #  load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               tau =  meta.res.indept$tau2,
                               I2 = meta.res.indept$I2,
                               pval.Q = meta.res.indept$pval.Q,
                               te.val.r = meta.res.indept$TE.random,
                               sete.val.r = meta.res.indept$seTE.random,
                               te.val.f = meta.res.indept$TE.fixed,
                               sete.val.f = meta.res.indept$seTE.fixed,
                               pval.r = meta.res.indept$pval.random,
                               pval.f = meta.res.indept$pval.fixed)

###############################################################################
## Five duplicate analysis
###############################################################################

id.res <- lapply(1:3, function(j){
  
  id.mod <- lapply(1:(length(study.id)-(j+3)), function(k){
    
    id1 <- study.id[k]
    id2 <- study.id[k+j]
    id3 <- study.id[k+j+1]
    id4 <- study.id[k+j+2]
    
    id.dup <-  sapply((k+(j+3)):length(study.id), function(i){
      
      cbind(id1, id2, id3, id4, study.id[i])
      
    })
    
    id.dup
    
  })
  
  id.mod
  
})

mod.id <- lapply(1:length(id.res), function(k){
  
  do.call(cbind, id.res[[k]])
  
})

id <- do.call(cbind, mod.id)

res.all <- lapply(1:ncol(id), function(m){
  
    id.dup <- id[,m]
    id.missing <- id.dup
    id.known <- study.id[!(study.id %in% id.missing)]
    mod.dat.exprs <- dat.exprs
   
    
    dat.exprs.dup <- lapply(1:length(id.known), function(k){
      
      mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
      
      mod.dat.exprs
      
    })
    
    lm.res <- lapply(1:length(dat.exprs.dup), function(i){
      
      dat.expr <- dat.exprs.dup[[i]]
      
      res <- lapply(1:length(study.id), function(k){
        
        y <- dat.aac[, study.id[k]]
        y <- scale(as.numeric(y), center = T, scale = T)[,1]
        x <- dat.expr[ ,study.id[k]]
        
        fit <- lm(y ~ as.numeric(x))
        coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
        r.squared <- summary(fit)$r.squared 
        df <- fit$df.residual
        coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
        rownames(coef.fit) <- NULL
        colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
        coef.fit
        
      })
      
      res <- do.call(rbind, res)
      df.id <- paste("n", res$df + 2, sep=" = ")
      rownames(res) <- paste(study.id, df.id, sep=", ")
      id.miss <- paste(id.missing[1], id.missing[2], 
                       id.missing[3], id.missing[4], id.missing[5], sep="-")
      res <- data.frame(study.unknown = id.miss, 
                        study.dup = id.known[i], res)
      res
    
  })
  
 
    meta.res <- lapply(1:length(lm.res), function(j){
      
      x <-  lm.res[[j]]$estimate
      y <-  lm.res[[j]]$std
      dat.xy <- data.frame(x, y)
      
      rownames(dat.xy) <- rownames(lm.res)
      
      res <- metagen(x, y, studlab = rownames(dat.xy), 
                     method.tau = "DL")
      
      res
      
    }) 
 
    
    id.missing.dup <- lapply(1:length(lm.res), function(j){
      
      data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                 id.duplicate = unique(lm.res[[j]]$study.dup))
      
    })
    
    id.miss.dup <- do.call(rbind, id.missing.dup)
    
    meta.result <- lapply(1:length(meta.res), function(k){
      
      data.frame(id.miss.dup[k,],
                 tau =  meta.res[[k]]$tau2,
                 I2 = meta.res[[k]]$I2,
                 pval.Q = meta.res[[k]]$pval.Q,
                 te.val.r = meta.res[[k]]$TE.random,
                 sete.val.r = meta.res[[k]]$seTE.random,
                 te.val.f = meta.res[[k]]$TE.fixed,
                 sete.val.f = meta.res[[k]]$seTE.fixed,
                 pval.r = meta.res[[k]]$pval.random,
                 pval.f = meta.res[[k]]$pval.fixed)
      
      
    })
    
    do.call(rbind, meta.result)
   
  })

meta.non.indept <- do.call(rbind, res.all)
dat.gene <- rbind( meta.indept,  meta.non.indept)
data.frame(gene_name = id.gene[t],
           drug = drugs[l], 
           dat.gene)

})  

dat.meta <- do.call(rbind, res.gene.drug) 
dat.meta

})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/five.dup.RData")


#############################################################
#############################################################
## six duplication
#############################################################
#############################################################


dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
 #   load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               tau =  meta.res.indept$tau2,
                               I2 = meta.res.indept$I2,
                               pval.Q = meta.res.indept$pval.Q,
                               te.val.r = meta.res.indept$TE.random,
                               sete.val.r = meta.res.indept$seTE.random,
                               te.val.f = meta.res.indept$TE.fixed,
                               sete.val.f = meta.res.indept$seTE.fixed,
                               pval.r = meta.res.indept$pval.random,
                               pval.f = meta.res.indept$pval.fixed)
    
###############################################################################
## six duplicate analysis
###############################################################################

res.all <- lapply(1:length(study.id), function(m){
  
  id.known <- study.id[m]
  id.missing <- study.id[study.id != id.known]
  mod.dat.exprs <- dat.exprs
  
  dat.exprs.dup <- lapply(1:length(id.known), function(k){
    
    mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
    
    mod.dat.exprs
    
  })
  
  
  lm.res <- lapply(1:length(dat.exprs.dup), function(i){
    
    dat.expr <- dat.exprs.dup[[i]]
    
    res <- lapply(1:length(study.id), function(k){
      
      y <- dat.aac[, study.id[k]]
      y <- scale(as.numeric(y), center = T, scale = T)[,1]
      x <- dat.expr[ ,study.id[k]]
      
      fit <- lm(y ~ as.numeric(x))
      coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
      r.squared <- summary(fit)$r.squared 
      df <- fit$df.residual
      coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
      rownames(coef.fit) <- NULL
      colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
      coef.fit
      
    })
    
    res <- do.call(rbind, res)
    df.id <- paste("n", res$df + 2, sep=" = ")
    rownames(res) <- paste(study.id, df.id, sep=", ")
    id.miss <- paste(id.missing[1], id.missing[2], 
                     id.missing[3], id.missing[4],
                     id.missing[5], id.missing[6], sep="-")
    res <- data.frame(study.unknown = id.miss, 
                      study.dup = id.known[i], res)
    res
    
  })
  
  
  meta.res <- lapply(1:length(lm.res), function(j){
    
    x <-  lm.res[[j]]$estimate
    y <-  lm.res[[j]]$std
    dat.xy <- data.frame(x, y)
    
    rownames(dat.xy) <- rownames(lm.res)
    
    res <- metagen(x, y, studlab = rownames(dat.xy), 
                   method.tau = "DL")
    
    res
    
  })
  
  id.missing.dup <- lapply(1:length(lm.res), function(j){
    
    data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
               id.duplicate = unique(lm.res[[j]]$study.dup))
    
  })
  
  id.miss.dup <- do.call(rbind, id.missing.dup)
  
  meta.result <- lapply(1:length(meta.res), function(k){
    
    data.frame(id.miss.dup[k,],
               tau =  meta.res[[k]]$tau2,
               I2 = meta.res[[k]]$I2,
               pval.Q = meta.res[[k]]$pval.Q,
               te.val.r = meta.res[[k]]$TE.random,
               sete.val.r = meta.res[[k]]$seTE.random,
               te.val.f = meta.res[[k]]$TE.fixed,
               sete.val.f = meta.res[[k]]$seTE.fixed,
               pval.r = meta.res[[k]]$pval.random,
               pval.f = meta.res[[k]]$pval.fixed)
    
  })
  
  do.call(rbind, meta.result)
  
})


meta.non.indept <- do.call(rbind, res.all)
dat.gene <- rbind( meta.indept,  meta.non.indept)
data.frame(gene_name = id.gene[t],
           drug = drugs[l], 
           dat.gene)

  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/six.dup.RData")


################################################################################################
################################################################################################
#################################### non-independent analyses: Bayesian ########################
################################################################################################
################################################################################################

load("/data/breast/DrugMeta.breast.zfixsample.RData")
load("/data/breast/DrugMeta.breast.fixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

id.gene <- c("CHMP7", "FBXO3", "DUOX2", "EGFR", "CD63", "ZNF143",
             "PSMB3", "C3orf52", "PITX3", "ERBB2", "CD63", "CD84",
             "IGF2BP3", "TGM3", "ABCG4", "S100A1", "COL11A1", "SCRT1")

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               I2 = meta.res.indept$I2,
                               te.val.r = meta.res.indept$TE.random,
                               ci.low = meta.res.indept$lower.random,
                               ci.up = meta.res.indept$upper.random)
    
    ##############################
    ## One duplicate analysis
    ##############################
    
    res.all <- lapply(1:length(study.id), function(m){
      
      id.missing <- study.id[m]
      id.known <- study.id[study.id != id.missing]
      mod.dat.exprs <- dat.exprs
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        mod.dat.exprs
        
      })
      
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        res <- data.frame(study.unknown = id.missing, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- bayesmeta(y=dat.xy$x, 
                         sigma=dat.xy$y,
                         mu.prior = c("mean"=NA, "sd"=NA), 
                         tau.prior="Jeffreys")
        
        res
        
      })
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   #tau =  meta.res[[k]]$summary["mean", "tau"], 
                   I2 = meta.res[[k]]$I2(tau = meta.res[[k]]$summary["mean", "tau"]),
                   te.val.r = meta.res[[k]]$summary["mean", "mu"],
                   ci.low = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[1],
                   ci.up = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[2])
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/one.dup.bayes.RData")

################################################################
################################################################
## two duplication
################################################################
#################################################################


dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               I2 = meta.res.indept$I2,
                               te.val.r = meta.res.indept$TE.random,
                               ci.low = meta.res.indept$lower.random,
                               ci.up = meta.res.indept$upper.random)
    
    ####################################################################
    ## two duplicate analysis
    ####################################################################
    
    id <- lapply(1:(length(study.id)-1), function(k){
      
      id1 <- study.id[k]
      
      
      id.dup <-  sapply((k+1):length(study.id), function(i){
        
        cbind(id1, study.id[i])
        
      })
      
      id.dup
      
    })
    
    
    res.all <- lapply(1:length(id), function(m){
      
      id.dup <- id[[m]]
      
      sub.res <-  lapply(1:ncol(id.dup), function(s){
        
        id.missing <- id.dup[, s]
        id.known <- study.id[!(study.id %in% id.missing)]
        mod.dat.exprs <- dat.exprs
        
        
        dat.exprs.dup <- lapply(1:length(id.known), function(k){
          
          mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
          
          mod.dat.exprs
          
        })
        
        lm.res <- lapply(1:length(dat.exprs.dup), function(i){
          
          dat.expr <- dat.exprs.dup[[i]]
          
          res <- lapply(1:length(study.id), function(k){
            
            y <- dat.aac[, study.id[k]]
            y <- scale(as.numeric(y), center = T, scale = T)[,1]
            x <- dat.expr[ ,study.id[k]]
            
            fit <- lm(y ~ as.numeric(x))
            coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
            r.squared <- summary(fit)$r.squared 
            df <- fit$df.residual
            coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
            rownames(coef.fit) <- NULL
            colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
            coef.fit
            
          })
          
          res <- do.call(rbind, res)
          df.id <- paste("n", res$df + 2, sep=" = ")
          rownames(res) <- paste(study.id, df.id, sep=", ")
          id.miss <- paste(id.missing[1], id.missing[2], sep="-")
          res <- data.frame(study.unknown = id.miss, 
                            study.dup = id.known[i], res)
          res
          
        })
        
        
        meta.res <- lapply(1:length(lm.res), function(j){
          
          x <-  lm.res[[j]]$estimate
          y <-  lm.res[[j]]$std
          dat.xy <- data.frame(x, y)
          
          rownames(dat.xy) <- rownames(lm.res)
          
          res <- bayesmeta(y=dat.xy$x, 
                           sigma=dat.xy$y,
                           mu.prior = c("mean"=NA, "sd"=NA), 
                           tau.prior="Jeffreys")
          
          res
          
        }) 
        
        
        id.missing.dup <- lapply(1:length(lm.res), function(j){
          
          data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                     id.duplicate = unique(lm.res[[j]]$study.dup))
          
        })
        
        id.miss.dup <- do.call(rbind, id.missing.dup)
        
        meta.result <- lapply(1:length(meta.res), function(k){
          
          data.frame(id.miss.dup[k,],
                     #tau =  meta.res[[k]]$summary["mean", "tau"], 
                     I2 = meta.res[[k]]$I2(tau = meta.res[[k]]$summary["mean", "tau"]),
                     te.val.r = meta.res[[k]]$summary["mean", "mu"],
                     ci.low = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[1],
                     ci.up = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[2])
          
        })
        
        do.call(rbind, meta.result)
        
      })
      
      
      
      do.call(rbind, sub.res)
      
      
    })
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/two.dup.bayes.RData")

#############################################################
#############################################################
## three duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               I2 = meta.res.indept$I2,
                               te.val.r = meta.res.indept$TE.random,
                               ci.low = meta.res.indept$lower.random,
                               ci.up = meta.res.indept$upper.random)
    
    ###################################################################
    ## Three duplicate analysis
    ####################################################################
    
    id.res <- lapply(1:5, function(j){
      
      id.mod <- lapply(1:(length(study.id)-(j+1)), function(k){
        
        id1 <- study.id[k]
        id2 <- study.id[k+j]
        
        id.dup <-  sapply((k+(j+1)):length(study.id), function(i){
          
          cbind(id1, id2, study.id[i])
          
        })
        
        id.dup
        
      })
      
      id.mod
      
    })
    
    mod.id <- lapply(1:length(id.res), function(k){
      
      do.call(cbind, id.res[[k]])
      
    })
    
    id <- do.call(cbind, mod.id)
    
    
    res.all <- lapply(1:ncol(id), function(m){
      
      id.dup <- id[,m]
      id.missing <- id.dup
      id.known <- study.id[!(study.id %in% id.missing)]
      mod.dat.exprs <- dat.exprs
      
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        
        mod.dat.exprs
        
      })
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        id.miss <- paste(id.missing[1], id.missing[2], id.missing[3], sep="-")
        res <- data.frame(study.unknown = id.miss, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- bayesmeta(y=dat.xy$x, 
                         sigma=dat.xy$y,
                         mu.prior = c("mean"=NA, "sd"=NA), 
                         tau.prior="Jeffreys")
        
        res
        
      }) 
      
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   #tau =  meta.res[[k]]$summary["mean", "tau"], 
                   I2 = meta.res[[k]]$I2(tau = meta.res[[k]]$summary["mean", "tau"]),
                   te.val.r = meta.res[[k]]$summary["mean", "mu"],
                   ci.low = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[1],
                   ci.up = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[2])
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/three.dup.bayes.RData")

#############################################################
#############################################################
## four duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               I2 = meta.res.indept$I2,
                               te.val.r = meta.res.indept$TE.random,
                               ci.low = meta.res.indept$lower.random,
                               ci.up = meta.res.indept$upper.random)
    
    ####################################################################
    ## Four duplicate analysis
    ####################################################################
    
    id.res <- lapply(1:4, function(j){
      
      id.mod <- lapply(1:(length(study.id)-(j+2)), function(k){
        
        id1 <- study.id[k]
        id2 <- study.id[k+j]
        id3 <- study.id[k+j+1]
        
        id.dup <-  sapply((k+(j+2)):length(study.id), function(i){
          
          cbind(id1, id2, id3, study.id[i])
          
        })
        
        id.dup
        
      })
      
      id.mod
      
    })
    
    mod.id <- lapply(1:length(id.res), function(k){
      
      do.call(cbind, id.res[[k]])
      
    })
    
    id <- do.call(cbind, mod.id)
    
    res.all <- lapply(1:ncol(id), function(m){
      
      id.dup <- id[,m]
      id.missing <- id.dup
      id.known <- study.id[!(study.id %in% id.missing)]
      mod.dat.exprs <- dat.exprs
      
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        
        mod.dat.exprs
        
      })
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        id.miss <- paste(id.missing[1], id.missing[2], 
                         id.missing[3], id.missing[4], sep="-")
        res <- data.frame(study.unknown = id.miss, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- bayesmeta(y=dat.xy$x, 
                         sigma=dat.xy$y,
                         mu.prior = c("mean"=NA, "sd"=NA), 
                         tau.prior="Jeffreys")
        
        res
        
      }) 
      
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   #tau =  meta.res[[k]]$summary["mean", "tau"], 
                   I2 = meta.res[[k]]$I2(tau = meta.res[[k]]$summary["mean", "tau"]),
                   te.val.r = meta.res[[k]]$summary["mean", "mu"],
                   ci.low = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[1],
                   ci.up = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[2])
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/four.dup.bayes.RData")

#############################################################
#############################################################
## five duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                               I2 = meta.res.indept$I2,
                               te.val.r = meta.res.indept$TE.random,
                               ci.low = meta.res.indept$lower.random,
                               ci.up = meta.res.indept$upper.random)
    
    ###############################################################################
    ## Five duplicate analysis
    ###############################################################################
    
    id.res <- lapply(1:3, function(j){
      
      id.mod <- lapply(1:(length(study.id)-(j+3)), function(k){
        
        id1 <- study.id[k]
        id2 <- study.id[k+j]
        id3 <- study.id[k+j+1]
        id4 <- study.id[k+j+2]
        
        id.dup <-  sapply((k+(j+3)):length(study.id), function(i){
          
          cbind(id1, id2, id3, id4, study.id[i])
          
        })
        
        id.dup
        
      })
      
      id.mod
      
    })
    
    mod.id <- lapply(1:length(id.res), function(k){
      
      do.call(cbind, id.res[[k]])
      
    })
    
    id <- do.call(cbind, mod.id)
    
    res.all <- lapply(1:ncol(id), function(m){
      
      id.dup <- id[,m]
      id.missing <- id.dup
      id.known <- study.id[!(study.id %in% id.missing)]
      mod.dat.exprs <- dat.exprs
      
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        
        mod.dat.exprs
        
      })
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        id.miss <- paste(id.missing[1], id.missing[2], 
                         id.missing[3], id.missing[4], id.missing[5], sep="-")
        res <- data.frame(study.unknown = id.miss, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- bayesmeta(y=dat.xy$x, 
                         sigma=dat.xy$y,
                         mu.prior = c("mean"=NA, "sd"=NA), 
                         tau.prior="Jeffreys")
        
        res
        
      }) 
      
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   #tau =  meta.res[[k]]$summary["mean", "tau"], 
                   I2 = meta.res[[k]]$I2(tau = meta.res[[k]]$summary["mean", "tau"]),
                   te.val.r = meta.res[[k]]$summary["mean", "mu"],
                   ci.low = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[1],
                   ci.up = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[2])
        
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/five.dup.bayes.RData")


#############################################################
#############################################################
## six duplication
#############################################################
#############################################################


dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.uhn  <- z.exprs.uhn.data[z.exprs.uhn.data$gene_name == id.gene[t], ][-1]
    sub.gray  <- z.exprs.gray.data[z.exprs.gray.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.uhn <- mod.aac.uhn[ ,colnames(mod.aac.uhn) == drugs[l] ]
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.gray <- mod.aac.gray[ ,colnames(mod.aac.gray) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.uhn <- scale(as.numeric(sub.uhn), center = T, scale = T)
    expr.gray <- scale(as.numeric(sub.gray), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            UHN = as.numeric(expr.uhn[,1]), 
                            GRAY = as.numeric(expr.gray[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          UHN = as.numeric( sub.aac.uhn), 
                          GRAY = as.numeric( sub.aac.gray), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
    ##############################
    ## Independent analysis
    #############################
    
    #load("/results/Meta/Breast/meta.breast.RData")
    
    meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                              meta.effect.size.res$drug == drugs[l] &
                                              meta.effect.size.res$method.meta == "DL", ]
    
    meta.indept <- data.frame(id.missing = "Indept", id.duplicate = "Indept",
                              I2 = meta.res.indept$I2,
                              te.val.r = meta.res.indept$TE.random,
                              ci.low = meta.res.indept$lower.random,
                              ci.up = meta.res.indept$upper.random)
    
    ###############################################################################
    ## Six duplicate analysis
    ###############################################################################
    
    res.all <- lapply(1:length(study.id), function(m){
      
      id.known <- study.id[m]
      id.missing <- study.id[study.id != id.known]
      mod.dat.exprs <- dat.exprs
      
      dat.exprs.dup <- lapply(1:length(id.known), function(k){
        
        mod.dat.exprs[, id.missing] <-  dat.exprs[, id.known[k]]
        
        mod.dat.exprs
        
      })
      
      
      lm.res <- lapply(1:length(dat.exprs.dup), function(i){
        
        dat.expr <- dat.exprs.dup[[i]]
        
        res <- lapply(1:length(study.id), function(k){
          
          y <- dat.aac[, study.id[k]]
          y <- scale(as.numeric(y), center = T, scale = T)[,1]
          x <- dat.expr[ ,study.id[k]]
          
          fit <- lm(y ~ as.numeric(x))
          coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
          r.squared <- summary(fit)$r.squared 
          df <- fit$df.residual
          coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
          rownames(coef.fit) <- NULL
          colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
          coef.fit
          
        })
        
        res <- do.call(rbind, res)
        df.id <- paste("n", res$df + 2, sep=" = ")
        rownames(res) <- paste(study.id, df.id, sep=", ")
        id.miss <- paste(id.missing[1], id.missing[2], 
                         id.missing[3], id.missing[4],
                         id.missing[5], id.missing[6], sep="-")
        res <- data.frame(study.unknown = id.miss, 
                          study.dup = id.known[i], res)
        res
        
      })
      
      
      meta.res <- lapply(1:length(lm.res), function(j){
        
        x <-  lm.res[[j]]$estimate
        y <-  lm.res[[j]]$std
        dat.xy <- data.frame(x, y)
        
        rownames(dat.xy) <- rownames(lm.res)
        
        res <- bayesmeta(y=dat.xy$x, 
                         sigma=dat.xy$y,
                         mu.prior = c("mean"=NA, "sd"=NA), 
                         tau.prior="Jeffreys")
        
        res
        
      })
      
      id.missing.dup <- lapply(1:length(lm.res), function(j){
        
        data.frame(id.missing = unique(lm.res[[j]]$study.unknown),
                   id.duplicate = unique(lm.res[[j]]$study.dup))
        
      })
      
      id.miss.dup <- do.call(rbind, id.missing.dup)
      
      meta.result <- lapply(1:length(meta.res), function(k){
        
        data.frame(id.miss.dup[k,],
                   #tau =  meta.res[[k]]$summary["mean", "tau"], 
                   I2 = meta.res[[k]]$I2(tau = meta.res[[k]]$summary["mean", "tau"]),
                   te.val.r = meta.res[[k]]$summary["mean", "mu"],
                   ci.low = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[1],
                   ci.up = meta.res[[k]]$post.interval(mu.level=0.95, method="central")[2])
        
      })
      
      do.call(rbind, meta.result)
      
    })
    
    
    meta.non.indept <- do.call(rbind, res.all)
    dat.gene <- rbind( meta.indept,  meta.non.indept)
    data.frame(gene_name = id.gene[t],
               drug = drugs[l], 
               dat.gene)
    
  })  
  
  dat.meta <- do.call(rbind, res.gene.drug) 
  dat.meta
  
})

dat.meta.dup <- do.call(rbind, dat.meta.dup)

save(dat.meta.dup, 
     file="/results/MetaNoIndept/Breast/six.dup.bayes.RData")




################################################################################################################
################################################################################################################
########################################### visualization plot #################################################
################################################################################################################
################################################################################################################

####################################################################
## Breast : duplication
####################################################################

load("/results/MetaNoIndept/Breast/one.dup.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Breast/two.dup.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Breast/three.dup.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Breast/four.dup.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

load("/results/MetaNoIndept/Breast/five.dup.RData")
dat5 <- data.frame(missing.GEx = "five", dat.meta.dup)

load("/results/MetaNoIndept/Breast/six.dup.RData")
dat6 <- data.frame(missing.GEx = "six", dat.meta.dup)

load("/results/Cor/cor.breast.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id <- unique(dat1$gene_name)

res.gene.drug.all <- lapply(1:length(drugs), function(l){
  print(l)
  res.gene.drug <- lapply(1:length(id), function(k){
    print(k)
    sub.dat1 <- dat1[dat1$gene_name == id[k] & dat1$drug == drugs[l], ]
    sub.dat1 <- data.frame(samples = nrow(sub.dat1)-1, sub.dat1)
    
    sub.dat2 <- dat2[dat2$gene_name == id[k] & dat2$drug == drugs[l], ]
    sub.dat2 <- data.frame(samples = nrow(sub.dat2)-1, sub.dat2)
    
    sub.dat3 <- dat3[dat3$gene_name == id[k] & dat3$drug == drugs[l], ]
    sub.dat3 <- data.frame(samples = nrow(sub.dat3)-1, sub.dat3)
    
    sub.dat4 <- dat4[dat4$gene_name == id[k]  & dat4$drug == drugs[l], ]
    sub.dat4 <- data.frame(samples = nrow(sub.dat4)-1, sub.dat4)
    
    sub.dat5 <- dat5[dat5$gene_name == id[k] & dat5$drug == drugs[l], ]
    sub.dat5 <- data.frame(samples = nrow(sub.dat5)-1, sub.dat5)
    
    sub.dat6 <- dat6[dat6$gene_name == id[k] & dat6$drug == drugs[l], ]
    sub.dat6 <- data.frame(samples = nrow(sub.dat6)-1, sub.dat6)
    
    sub.data <- rbind(sub.dat1, sub.dat2, sub.dat3,
                      sub.dat4, sub.dat5, sub.dat6)
    
    indept.data <- sub.data[sub.data$id.duplicate == "Indept", ]
    
    sub.data$bias.te.f <- sub.data$te.val.f - indept.data$te.val.f[1]
    sub.data$bias.te.r <- sub.data$te.val.r - indept.data$te.val.r[1]
    
    data <- sub.data[sub.data$id.duplicate != "Indept", ]
    sub.data$id.duplicate <- factor(sub.data$id.duplicate)
    id.miss <- unique(data$missing.GEx)
    id.samples <- unique(data$samples)
    
    dat <- data[data$gene == id[k], ]
    
    mean.dat <- lapply(1:length(id.miss), function(j){
      
      bias.te.f.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.f))
      sd.te.f.val <- sd(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.f))
      
      bias.te.r.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))
      sd.te.r.val <- sd(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))  
      
      tau.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$tau)
      I2.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$I2)
      
      pval.Q.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$pval.Q)
      
      data.frame(id = id.miss[j], samples = id.samples[j], 
                 drug = unique(sub.data$drug), gene = id[k], 
                 I2.indept = unique(indept.data$I2),
                 pval.Q = unique(indept.data$pval.Q),
                 tau.mean = tau.val,
                 I2.mean = I2.val,
                 pval.Q.mean = pval.Q.val,
                 bias.te.f.val, bias.te.r.val,
                 sd.te.f.val, sd.te.r.val)
      
    })
    
    
    res.all <- do.call(rbind, mean.dat)
    res.all
    
  })
  
  do.call(rbind, res.gene.drug)
  
})

res.gene.drug.all <- do.call(rbind, res.gene.drug.all)


cor.bias.res.all <- lapply(1:length(drugs), function(l){
  
  cor.bias.res <- lapply(1:length(id), function(k){
    
    sub.cor.res <- cor.res[cor.res$gene == id[k], ]
    sub.res <- res.gene.drug.all[res.gene.drug.all$gene == id[k] & 
                                   res.gene.drug.all$drug == drugs[l], ]
    I2.mean <- mean(sub.res$I2.mean)
    bias.mean <- mean(sub.res$bias.te.r.val)
    data.frame(gene = unique(sub.res$gene), 
               drug = unique(sub.res$drug), 
               bias.mean = bias.mean, 
               I2.mean = I2.mean, 
               I2 = sub.res$I2.indept[1],
               sub.cor.res[, 2:3])
    
  })
  
  do.call(rbind, cor.bias.res)
  
})

cor.bias.res.all <- do.call(rbind, cor.bias.res.all)

save(cor.bias.res.all,  res.gene.drug.all ,
     file="/results/MetaNoIndept/Breast/cor.bias.res.RData")

#########################################################
##  Breast: Top genes to plot
#########################################################

load("/results/MetaNoIndept/Breast/cor.bias.res.RData")
#load("/results/Meta/Breast/meta.breast.RData")

drugs <- "Erlotinib"
id <- c("CHMP7", "FBXO3", "DUOX2", "EGFR", "CD63", "ZNF143", 
        "PSMB3",  "C3orf52", "PITX3", "ERBB2", "CD84", "IGF2BP3", 
        "TGM3", "ABCG4", "S100A1", "COL11A1", "SCRT1")

gene.id <- id[1:6]
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & meta.effect.size.res$method.meta == "DL", ]
res <- res[order(res$gene.name), ]

cor.res <- cor.bias.res.all[cor.bias.res.all$drug == drugs, ]
cor.res <- cor.res[order(cor.res$gene), ]
cor.res <- cor.res[cor.res$gene %in% gene.id, ]

res <- res[res$gene.name %in% cor.res$gene, ]

dat <- cbind(res, cor.res[, -c(1:2)])

dat$label <- NA
dat$label <- sapply(1:nrow(dat), function(k){
  
  cor.val <- dat$median.cor[k]
  if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
    if(abs(cor.val) > 0.7){  label.val <- "High"}else{
      label.val <- "Medium"
    }
  }
  
  label.val
  
})

dat$label <- factor(dat$label)
sig <- dat[dat$padj < 0.05, ]

sig.high.i2 <- sig[sig$I2 > 0.50 & sig$pval.Q < 0.1, ]
dat.high.i2 <- sig.high.i2[order(sig.high.i2$median.cor, decreasing = TRUE), ]
dat.high.i2 <- data.frame(dat.high.i2, I2.id = "I2 > 50%")

sig.low.i2 <- sig[!(sig$I2 > 0.50 & sig$pval.Q < 0.1), ]
dat.low.i2 <- sig.low.i2[order(sig.low.i2$median.cor, decreasing = TRUE), ]
dat.low.i2 <- rbind(dat[dat$gene.name == "EGFR", ],
                    dat.low.i2)  
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug ==drugs, ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(3,6,5,2,1,4)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(2,5,3,4,1,6)])

p <-ggplot(data.plot, 
           aes(y = as.numeric(bias.te.r.val), x = id, fill= I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  ylab("MAD pooled effect size, Erlotinib") +
  xlab("") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/Fig5_Erlotinib.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()


drugs <- "Lapatinib"
gene.id <- c(id[7:11], "CD63")
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & meta.effect.size.res$method.meta == "DL", ]
res <- res[order(res$gene.name), ]

cor.res <- cor.bias.res.all[cor.bias.res.all$drug == drugs, ]
cor.res <- cor.res[order(cor.res$gene), ]
cor.res <- cor.res[cor.res$gene %in% gene.id, ]

res <- res[res$gene.name %in% cor.res$gene, ]

dat <- cbind(res, cor.res[, -c(1:2)])

dat$label <- NA
dat$label <- sapply(1:nrow(dat), function(k){
  
  cor.val <- dat$median.cor[k]
  if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
    if(abs(cor.val) > 0.7){  label.val <- "High"}else{
      label.val <- "Medium"
    }
  }
  
  label.val
  
})

dat$label <- factor(dat$label)
dat$label <- factor(dat$label, levels = levels(dat$label)[c(1,3,2)])
sig <- dat[dat$padj < 0.05, ]

sig.high.i2 <- sig[sig$I2 > 0.50 & sig$pval.Q < 0.1, ]
dat.high.i2 <- sig.high.i2[order(sig.high.i2$median.cor, decreasing = TRUE), ]
dat.high.i2 <- data.frame(dat.high.i2, I2.id = "I2 > 50%")

sig.low.i2 <- sig[!(sig$I2 > 0.50 & sig$pval.Q < 0.1), ]
dat.low.i2 <- sig.low.i2[order(sig.low.i2$median.cor, decreasing=TRUE), ]
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug ==drugs, ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(3,6,5,2,1,4)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(6,1,5,4,2,3)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  ylab("MAD pooled effect size, Lapatinib") +
  xlab("") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/Fig5_Lapatinib.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)


print(p)

dev.off()

drugs <- "Paclitaxel"
gene.id <- id[12:17]
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & meta.effect.size.res$method.meta == "DL", ]
res <- res[order(res$gene.name), ]

cor.res <- cor.bias.res.all[cor.bias.res.all$drug == drugs, ]
cor.res <- cor.res[order(cor.res$gene), ]
cor.res <- cor.res[cor.res$gene %in% gene.id, ]

res <- res[res$gene.name %in% cor.res$gene, ]

dat <- cbind(res, cor.res[, -c(1:2)])

dat$label <- NA
dat$label <- sapply(1:nrow(dat), function(k){
  
  cor.val <- dat$median.cor[k]
  if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
    if(abs(cor.val) > 0.7){  label.val <- "High"}else{
      label.val <- "Medium"
    }
  }
  
  label.val
  
})

dat$label <- factor(dat$label)
dat$label <- factor(dat$label, levels = levels(dat$label)[c(1,3,2)])
sig <- dat[dat$padj < 0.05, ]

sig.high.i2 <- sig[sig$I2 > 0.50 & sig$pval.Q < 0.1, ]
dat.high.i2 <- sig.high.i2[order(sig.high.i2$median.cor, decreasing = TRUE), ]
dat.high.i2 <- data.frame(dat.high.i2, I2.id = "I2 > 50%")

sig.low.i2 <- sig[!(sig$I2 > 0.50 & sig$pval.Q < 0.1), ]
dat.low.i2 <- sig.low.i2[order(sig.low.i2$median.cor, decreasing = TRUE), ]
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug == drugs , ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(3,6,5,2,1,4)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(3,6,1,4,2,5)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  ylab("MAD pooled effect size, Paclitaxel") +
  xlab("") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/Fig5_Paclitaxel.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

###################################################
##  Breast: Correlation vs bias
##################################################

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
dat <- lapply(1:length(drugs), function(l){
  
  sub.dat <- cor.bias.res.all[cor.bias.res.all$drug == drugs[l], ]
  sub.dat$label <- NA
  sub.dat$label <- sapply(1:nrow(sub.dat), function(k){
    
    cor.val <- sub.dat$median.cor[k]
    if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
      if(abs(cor.val) > 0.7){  label.val <- "High"}else{
        label.val <- "Medium"
      }
    }
    
    label.val
    
  })
  
  sub.dat$label <- factor(sub.dat$label)
  sub.dat$label <- factor(sub.dat$label, levels = levels(sub.dat$label)[c(2,3,1)])
  sub.dat
  
  
})

dat <- do.call(rbind, dat)

p <-ggplot(dat, aes(label, bias.mean, fill = label)) +
  geom_violin(width = 0.8) +
  geom_boxplot(width=0.08) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  scale_fill_manual(values=c("indianred2", "cadetblue3", "#E69F00")) + 
  facet_wrap( ~ drug, ncol=3) +
  xlab("") +
  ylab("average of MAD of pooled effect size") +
  theme(axis.text.x=element_text(size=14,  face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=14),
        strip.text= element_text(size=12, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE)

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/Fig5_cor_bias.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)
print(p)

dev.off()

################################################################
## Breast: Jaccard Index
################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

# NOTE: the following code was commented due to taking time to get the non-independent results across all genes on CodeOcean. 
# We only run for subset of genes on CodeOcean. We already have the results across all genes on HPC.

#load("/results/Meta/Breast/meta.breast.RData")

#load("/results/MetaNoIndept/Breast/one.dup.RData")
#dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

#load("/results/MetaNoIndept/Breast/two.dup.RData")
#dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

#load("/results/MetaNoIndept/Breast/three.dup.RData")
#dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

#load("/results/MetaNoIndept/Breast/four.dup.RData")
#dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

#load("/results/MetaNoIndept/Breast/five.dup.RData")
#dat5 <- data.frame(missing.GEx = "five", dat.meta.dup)

#load("/results/MetaNoIndept/Breast/six.dup.RData")
#dat6 <- data.frame(missing.GEx = "six", dat.meta.dup)

#dat.file <- rbind(dat1, dat2, dat3, dat4, dat5, dat6) 
#drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

#res.JC <- lapply(1:length(drugs), function(l){

#  dat.indept <- meta.effect.size.res[meta.effect.size.res$drug == drugs[l]
#                                     & meta.effect.size.res$method.meta == "DL", ]

#  dat.indept <- dat.indept[order(dat.indept$pval.random), ]
#  top.sig.gene <- dat.indept$gene.name[1:100]

#  data <- dat.file[dat.file$drug == drugs[l], ]
#  id.missing.GEx <- unique(data$missing.GEx)

#  sub.data <-  lapply(1:length(id.missing.GEx), function(i){

#    dat <- data[data$missing.GEx == id.missing.GEx[i], ]

#    id.missing <- unique(dat$id.missing)[-1]

#    mod.dat <- lapply(1:length(id.missing), function(k){

#      sub.dat <- dat[dat$id.missing == id.missing[k], ]
#      id.duplicate <- unique(sub.dat$id.duplicate)

#      dat.dup <- lapply(1:length(id.duplicate), function(j){

#        sub.dat.dup <- sub.dat[sub.dat$id.duplicate == id.duplicate[j], ]
#        sub.dat.dup[order(sub.dat.dup$pval.r), ][1:100, ]

#      })

#      do.call(rbind, dat.dup)

#    })

#    do.call(rbind, mod.dat)

#  })

#  sub.dat <- do.call(rbind, sub.data)


#  JI.similarity.res <- lapply(1:length(id.missing.GEx), function(i){

#    id.missing <- unique(sub.dat[sub.dat$missing.GEx == id.missing.GEx[i], "id.missing"])

#    similarity.res <- lapply(1:length(id.missing), function(k){

#      id <- unique(sub.dat[sub.dat$id.missing == id.missing[k] &
#                             sub.dat$missing.GEx == id.missing.GEx[i], "id.duplicate"])

#      index.val <- sapply(1:length(id), function(j){

#        jaccard(sub.dat[sub.dat$missing.GEx == id.missing.GEx[i] &
#                          sub.dat$id.duplicate == id[j] & 
#                          sub.dat$id.missing == id.missing[k], "gene_name"],
#                dat.indept$gene.name[1:100])

#      }) 

#      data.frame(missing.GEx = id.missing.GEx[i],
#                 drug = sub.dat$drug[1],
#                 id.missing = id.missing[k],
#                 id.duplicate = id, 
#                 JI = index.val)

#    })

#    do.call(rbind,  similarity.res)

#  })

#  do.call(rbind, JI.similarity.res)


# })

# dat.plot <- do.call(rbind, res.JC)
# dat.plot$missing.GEx <- factor(dat.plot$missing.GEx)
# dat.plot$missing.GEx <- factor(dat.plot$missing.GEx,
#                               levels = levels(dat.plot$missing.GEx)[c(3,6,5,2,1,4)])

#p <-ggplot(dat.plot , aes(x= missing.GEx, JI, fill = missing.GEx)) +
#  geom_violin(width = 0.8) +
#  geom_boxplot(width=0.06) +
#  stat_summary(fun.y=median, geom="point", size=1.5, color="black") +
#  scale_fill_brewer(palette="Blues")  +
#  facet_wrap( ~ drug, ncol=1) +
#  xlab("") +
#  ylab("Jaccard similarity index") +
#  theme(axis.text.x=element_text(size=14,  face="bold"),
#        axis.title=element_text(size=14, face="bold"),
#        axis.text.y=element_text(size=12),
#        strip.text= element_text(size=12, face="bold"),
#        panel.grid.major = element_blank(), 
#        panel.grid.minor = element_blank(),
#        panel.background = element_blank(),
#        plot.background = element_blank(), 
#        axis.line = element_line(colour = "black")) +
#  guides(fill=FALSE)


#Cairo::Cairo(
#  15, #length
#  15, #width
#  file = "/results/MetaNoIndept/Breast/Fig7A.jpeg",
#  type = "jpeg", 
#  bg = "transparent", 
#  dpi = 300,
#  units = "cm"
#)

#print(p)

#dev.off()

###############################################################
## Breast: Trend test
###############################################################

# NOTE: only run for selected genes on CodeOcean, while we did trend test across all genes on HPC.

load("/results/MetaNoIndept/Breast/cor.bias.res.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

trend.res <- lapply(1:length(drugs), function(l){
  
  dat <- res.gene.drug.all[res.gene.drug.all$drug == drugs[l], ]
  id <- unique(dat$gene)
  
  trend.pval <- lapply(1:length(id), function(k){
    
    sub.dat <- dat[dat$gene == id[k], ]
    trend.fit.kendal <- mk.test(sub.dat$bias.te.r.val, alternative = "greater")
    data.frame(gene = id[k], drug = drugs[l], pval.kendal = trend.fit.kendal$p.value)
    
  })
  
  do.call(rbind, trend.pval)
  
})

trend.res <- do.call(rbind,  trend.res)

write.csv(trend.res, file="/results/MetaNoIndept/Breast/trend.test.csv")

################################################################################################################
################################################################################################################
########################################### visualization plot: Bayesian #######################################
################################################################################################################
################################################################################################################

####################################################################
## Breast : duplication
####################################################################
load("/results/MetaNoIndept/Breast/one.dup.bayes.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Breast/two.dup.bayes.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Breast/three.dup.bayes.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Breast/four.dup.bayes.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

load("/results/MetaNoIndept/Breast/five.dup.bayes.RData")
dat5 <- data.frame(missing.GEx = "five", dat.meta.dup)

load("/results/MetaNoIndept/Breast/six.dup.bayes.RData")
dat6 <- data.frame(missing.GEx = "six", dat.meta.dup)

load("/results/Cor/cor.breast.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id <- unique(dat1$gene_name)

res.gene.drug.all <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id), function(k){
    
    sub.dat1 <- dat1[dat1$gene_name == id[k] & dat1$drug == drugs[l], ]
    sub.dat1 <- data.frame(samples = nrow(sub.dat1)-1, sub.dat1)
    
    sub.dat2 <- dat2[dat2$gene_name == id[k] & dat2$drug == drugs[l], ]
    sub.dat2 <- data.frame(samples = nrow(sub.dat2)-1, sub.dat2)
    
    sub.dat3 <- dat3[dat3$gene_name == id[k] & dat3$drug == drugs[l], ]
    sub.dat3 <- data.frame(samples = nrow(sub.dat3)-1, sub.dat3)
    
    sub.dat4 <- dat4[dat4$gene_name == id[k]  & dat4$drug == drugs[l], ]
    sub.dat4 <- data.frame(samples = nrow(sub.dat4)-1, sub.dat4)
    
    sub.dat5 <- dat5[dat5$gene_name == id[k] & dat5$drug == drugs[l], ]
    sub.dat5 <- data.frame(samples = nrow(sub.dat5)-1, sub.dat5)
    
    sub.dat6 <- dat6[dat6$gene_name == id[k] & dat6$drug == drugs[l], ]
    sub.dat6 <- data.frame(samples = nrow(sub.dat6)-1, sub.dat6)
    
    sub.data <- rbind(sub.dat1, sub.dat2, sub.dat3,
                      sub.dat4, sub.dat5, sub.dat6)
    
    indept.data <- sub.data[sub.data$id.duplicate == "Indept", ]
    indept.data$ci <- indept.data$ci.up - indept.data$ci.low
    
    sub.data$bias.te.r <- sub.data$te.val.r - indept.data$te.val.r[1]
    sub.data$ci <- sub.data$ci.up - sub.data$ci.low
    sub.data$ratio.ci <-  sub.data$ci / indept.data$ci[1]
    
    data <- sub.data[sub.data$id.duplicate != "Indept", ]
    sub.data$id.duplicate <- factor(sub.data$id.duplicate)
    id.miss <- unique(data$missing.GEx)
    id.samples <- unique(data$samples)
    
    dat <- data[data$gene == id[k], ]
    
    mean.dat <- lapply(1:length(id.miss), function(j){
      
      bias.te.r.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))
      I2.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$I2)
      
      ci.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$ratio.ci)
      
      data.frame(id = id.miss[j], samples = id.samples[j], 
                 drug = unique(sub.data$drug), gene = id[k], 
                 I2.indept = unique(indept.data$I2),
                 I2.mean = I2.val,
                 bias.te.r.val, 
                 ci.val)
      
    })
    
    
    res.all <- do.call(rbind, mean.dat)
    res.all
    
  })
  
  do.call(rbind, res.gene.drug)
  
})

res.gene.drug.all <- do.call(rbind, res.gene.drug.all)


cor.bias.res.all <- lapply(1:length(drugs), function(l){
  
  cor.bias.res <- lapply(1:length(id), function(k){
    
    sub.cor.res <- cor.res[cor.res$gene == id[k], ]
    sub.res <- res.gene.drug.all[res.gene.drug.all$gene == id[k] & 
                                   res.gene.drug.all$drug == drugs[l], ]
    I2.mean <- mean(sub.res$I2.mean)
    bias.mean <- mean(sub.res$bias.te.r.val)
    data.frame(gene = unique(sub.res$gene), 
               drug = unique(sub.res$drug), 
               bias.mean = bias.mean, 
               I2.mean = I2.mean, 
               I2 = sub.res$I2.indept[1],
               sub.cor.res[, 2:3])
    
  })
  
  do.call(rbind, cor.bias.res)
  
})

cor.bias.res.all <- do.call(rbind, cor.bias.res.all)

save(cor.bias.res.all,  res.gene.drug.all ,
     file="/results/MetaNoIndept/Breast/cor.bias.bayes.res.RData")

#########################################################
##  Breast: Top genes to plot
#########################################################

load("/results/MetaNoIndept/Breast/cor.bias.bayes.res.RData")
#load("/results/Meta/Breast/meta.breast.RData")

drugs <- "Erlotinib"
id <- c("CHMP7", "FBXO3", "DUOX2", "EGFR", "CD63",    "ZNF143", "PSMB3",
        "C3orf52", "PITX3", "ERBB2", "CD84", "IGF2BP3", "TGM3", 
        "ABCG4", "S100A1", "COL11A1", "SCRT1")

gene.id <- id[1:6]
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & meta.effect.size.res$method.meta == "DL", ]
res <- res[order(res$gene.name), ]

cor.res <- cor.bias.res.all[cor.bias.res.all$drug == drugs, ]
cor.res <- cor.res[order(cor.res$gene), ]
cor.res <- cor.res[cor.res$gene %in% gene.id, ]

res <- res[res$gene.name %in% cor.res$gene, ]

dat <- cbind(res, cor.res[, -c(1:2)])

dat$label <- NA
dat$label <- sapply(1:nrow(dat), function(k){
  
  cor.val <- dat$median.cor[k]
  if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
    if(abs(cor.val) > 0.7){  label.val <- "High"}else{
      label.val <- "Medium"
    }
  }
  
  label.val
  
})

dat$label <- factor(dat$label)
sig <- dat[dat$padj < 0.05, ]

sig.high.i2 <- sig[sig$I2 > 0.50 & sig$pval.Q < 0.1, ]
dat.high.i2 <- sig.high.i2[order(sig.high.i2$median.cor, decreasing = TRUE), ]
dat.high.i2 <- data.frame(dat.high.i2, I2.id = "I2 > 50%")

sig.low.i2 <- sig[!(sig$I2 > 0.50 & sig$pval.Q < 0.1), ]
dat.low.i2 <- sig.low.i2[order(sig.low.i2$median.cor, decreasing = TRUE), ]
dat.low.i2 <- rbind(dat[dat$gene.name == "EGFR", ],
                    dat.low.i2)  
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug ==drugs  , ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(3,6,5,2,1,4)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(2,5,3,4,1,6)])

p <-ggplot(data.plot, 
           aes(y = as.numeric(bias.te.r.val), x = id, fill= I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  # ggtitle(plot.id) +
  ylab("MAD pooled effect size, Erlotinib") +
  xlab("") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/SFig12_Erlotinib.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()


drugs <- "Lapatinib"
gene.id <- c(id[7:11], "CD63")
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & meta.effect.size.res$method.meta == "DL", ]
res <- res[order(res$gene.name), ]

cor.res <- cor.bias.res.all[cor.bias.res.all$drug == drugs, ]
cor.res <- cor.res[order(cor.res$gene), ]
cor.res <- cor.res[cor.res$gene %in% gene.id, ]

res <- res[res$gene.name %in% cor.res$gene, ]

dat <- cbind(res, cor.res[, -c(1:2)])

dat$label <- NA
dat$label <- sapply(1:nrow(dat), function(k){
  
  cor.val <- dat$median.cor[k]
  if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
    if(abs(cor.val) > 0.7){  label.val <- "High"}else{
      label.val <- "Medium"
    }
  }
  
  label.val
  
})

dat$label <- factor(dat$label)
dat$label <- factor(dat$label, levels = levels(dat$label)[c(1,3,2)])
sig <- dat[dat$padj < 0.05, ]

sig.high.i2 <- sig[sig$I2 > 0.50 & sig$pval.Q < 0.1, ]
dat.high.i2 <- sig.high.i2[order(sig.high.i2$median.cor, decreasing = TRUE), ]
dat.high.i2 <- data.frame(dat.high.i2, I2.id = "I2 > 50%")

sig.low.i2 <- sig[!(sig$I2 > 0.50 & sig$pval.Q < 0.1), ]
dat.low.i2 <- sig.low.i2[order(sig.low.i2$median.cor, decreasing=TRUE), ]
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug ==drugs , ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(3,6,5,2,1,4)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(6,1,5,4,2,3)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  # ggtitle(plot.id) +
  ylab("MAD pooled effect size, Lapatinib") +
  xlab("") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())


Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/SFig12_Lapatinib.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

drugs <- "Paclitaxel"
gene.id <- id[12:17]
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & meta.effect.size.res$method.meta == "DL", ]
res <- res[order(res$gene.name), ]

cor.res <- cor.bias.res.all[cor.bias.res.all$drug == drugs, ]
cor.res <- cor.res[order(cor.res$gene), ]
cor.res <- cor.res[cor.res$gene %in% gene.id, ]

res <- res[res$gene.name %in% cor.res$gene, ]

dat <- cbind(res, cor.res[, -c(1:2)])

dat$label <- NA
dat$label <- sapply(1:nrow(dat), function(k){
  
  cor.val <- dat$median.cor[k]
  if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
    if(abs(cor.val) > 0.7){  label.val <- "High"}else{
      label.val <- "Medium"
    }
  }
  
  label.val
  
})

dat$label <- factor(dat$label)
dat$label <- factor(dat$label, levels = levels(dat$label)[c(1,3,2)])
sig <- dat[dat$padj < 0.05, ]

sig.high.i2 <- sig[sig$I2 > 0.50 & sig$pval.Q < 0.1, ]
dat.high.i2 <- sig.high.i2[order(sig.high.i2$median.cor, decreasing = TRUE), ]
dat.high.i2 <- data.frame(dat.high.i2, I2.id = "I2 > 50%")

sig.low.i2 <- sig[!(sig$I2 > 0.50 & sig$pval.Q < 0.1), ]
dat.low.i2 <- sig.low.i2[order(sig.low.i2$median.cor, decreasing = TRUE), ]
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug ==drugs, ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(3,6,5,2,1,4)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(3,6,1,4,2,5)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  # ggtitle(plot.id) +
  ylab("MAD pooled effect size, Paclitaxel") +
  xlab("") +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=10, face = "bold"),
        strip.text = element_text(size=10, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())


Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaNoIndept/Breast/SFig12_Paclitaxel.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)


print(p)

dev.off()


###################################################
##  Breast: 95% CI DL vs Bayes credible interval
##################################################

load("/results/MetaNoIndept/Breast/cor.bias.bayes.res.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

data.plot <- lapply(1:length(drugs), function(l){
  
  
  bayes.one <- mean(res.gene.drug.all[res.gene.drug.all$id == "one" & res.gene.drug.all$drug == drugs[l],
                                      "ci.val"])
  bayes.two <- mean(res.gene.drug.all[res.gene.drug.all$id == "two" & res.gene.drug.all$drug == drugs[l], 
                                      "ci.val"])
  bayes.three <- mean(res.gene.drug.all[res.gene.drug.all$id == "three" & res.gene.drug.all$drug == drugs[l],
                                        "ci.val"])
  bayes.four <- mean(res.gene.drug.all[res.gene.drug.all$id == "four" & res.gene.drug.all$drug == drugs[l], 
                                       "ci.val"])
  bayes.five <- mean(res.gene.drug.all[res.gene.drug.all$id == "five" & res.gene.drug.all$drug == drugs[l], 
                                       "ci.val"])
  bayes.six <- mean(res.gene.drug.all[res.gene.drug.all$id == "six" & res.gene.drug.all$drug == drugs[l], 
                                      "ci.val"])
  
  load("/results/MetaNoIndept/Breast/one.dup.bayes.RData")
  dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)
  dat1 <- dat1[dat1$id.missing != "Indept", ]
  dat1 <- dat1[dat1$drug == drugs[l], ]
  dl.one <- mean(dat1$ci.up - dat1$ci.low)
  
  load("/results/MetaNoIndept/Breast/two.dup.bayes.RData")
  dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)
  dat2 <- dat2[dat2$id.missing != "Indept", ]
  dat2 <- dat2[dat2$drug == drugs[l], ]
  dl.two <- mean(dat2$ci.up - dat2$ci.low)
  
  load("/results/MetaNoIndept/Breast/three.dup.bayes.RData")
  dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)
  dat3 <- dat3[dat3$id.missing != "Indept", ]
  dat3 <- dat3[dat3$drug == drugs[l], ]
  dl.three <- mean(dat3$ci.up - dat3$ci.low)
  
  load("/results/MetaNoIndept/Breast/four.dup.bayes.RData")
  dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)
  dat4 <- dat4[dat4$id.missing != "Indept", ]
  dat4 <- dat4[dat4$drug == drugs[l], ]
  dl.four <- mean(dat4$ci.up - dat4$ci.low)
  
  load("/results/MetaNoIndept/Breast/five.dup.bayes.RData")
  dat5 <- data.frame(missing.GEx = "five", dat.meta.dup)
  dat5 <- dat5[dat5$id.missing != "Indept", ]
  dat5 <- dat5[dat5$drug == drugs[l], ]
  dl.five <- mean(dat5$ci.up - dat5$ci.low)
  
  load("/results/MetaNoIndept/Breast/six.dup.bayes.RData")
  dat6 <- data.frame(missing.GEx = "six", dat.meta.dup)
  dat6 <- dat6[dat6$id.missing != "Indept", ]
  dat6 <- dat6[dat6$drug == drugs[l], ]
  dl.six <- mean(dat6$ci.up - dat6$ci.low)
  
  data.frame(y = c(dl.one, dl.two, dl.three, dl.four, dl.five, dl.six,
                   bayes.one, bayes.two, bayes.three, bayes.four, bayes.five, bayes.six),
             drug = drugs[l],
             group = rep(c("one", "two", "three", "four", "five", "six"), 2), 
             method = c(rep("DerSimonian-Laird (DL)", 6), rep("Jeffreys", 6)))
})

dat <- do.call(rbind, data.plot)
dat$group <- factor(dat$group)
dat$group <- factor(dat$group, levels = levels(dat$group)[c(3,6,5,2,1,4)])

p <-ggplot(dat , aes(y =  y, x = group, fill = method) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  scale_fill_manual(values=c("maroon", "seashell3")) +
  facet_wrap(. ~ drug) +
  ylab("average of length of 95% interval") +
  xlab("") +
  theme(axis.text.x=element_text(size=12,  face="bold"),
        axis.title=element_text(size=14,face="bold"),
        axis.text.y=element_text(size=12, face = "bold"),
        strip.text = element_text(size=12, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="bottom",
        legend.text = element_text(size = 12, face="bold"),
        legend.title = element_blank())


Cairo::Cairo(
  20, #length
  14, #width
  file = "/results/MetaNoIndept/Breast/SFig12_CI_DL_Bayes.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()
