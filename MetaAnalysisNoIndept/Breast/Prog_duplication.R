source("/code/MetaAnalysis/Prog_MetaAnalysis.run.R")

#dir.create("/results/Meta")
dir.create("/results/MetaNoIndept")
dir.create("/results/MetaNoIndept/Breast")
################################################

load("/data/breast/DrugMeta.breast.zfixsample.RData")
load("/data/breast/DrugMeta.breast.fixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

id.gene <- z.exprs.gcsi.data$gene_name

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
    
    load("/results/Meta/Breast/meta.breast.RData")
    
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
    
     load("/results/Meta/Breast/meta.breast.RData")
    
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
    
    load("/results/Meta/Breast/meta.breast.RData")
    
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
    
    load("/results/Meta/Breast/meta.breast.RData")
    
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
    
    load("/results/Meta/Breast/meta.breast.RData")
    
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
    
    load("/results/Meta/Breast/meta.breast.RData")
    
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
## One duplicate analysis
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


