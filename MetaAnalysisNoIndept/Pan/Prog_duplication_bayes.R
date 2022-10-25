source("/code/MetaAnalysis/Prog_MetaAnalysis_fun.R")

#dir.create("/results/Meta")
dir.create("/results/MetaNoIndept")
dir.create("/results/MetaNoIndept/Pan")

################################################

load("/data/pan/DrugMeta.pan.zfixsample.RData")
load("/data/pan/DrugMeta.pan.fixsample.mod.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

id.gene <- c("C1orf116", "SPRR1B", "KRT1", "EGFR", "COX7A2", "ROS1" ,
        "TDRD1", "ZNF365", "KRT1", "ERBB2", "MAPK7", "CYP2A13", 
        "C4orf19", "ZNF688", "CPSF6", "S100A1", "NOC3L", "CD244")

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
####################################################################
## design of study
####################################################################
ccl.dat <- pan.ccl.ccle
ccl.dat$tissueid <- factor(ccl.dat$tissueid)
group.id <- unique(ccl.dat$tissueid)

group.tissue <- lapply(1:length(group.id), function(k){
  
  ifelse(ccl.dat$tissueid == group.id[k], 1, 0)
  
})

group.tissue <- do.call(cbind, group.tissue)
colnames( group.tissue) <- group.id

####################################################################
## Independent analysis
####################################################################

load("/results/Meta/Pan/meta.pan.RData")

meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                          meta.effect.size.res$drug == drugs[l] &
                                          meta.effect.size.res$method.meta == "DL", ]

meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                           I2 = meta.res.indept$I2,
                           te.val.r = meta.res.indept$TE.random,
                           ci.low = meta.res.indept$lower.random,
                           ci.up = meta.res.indept$upper.random)

res.all <- lapply(1:length(study.id), function(m){
  
  id.missing <- study.id[m]
  id.known <- study.id[study.id != id.missing]
  mod.dat.exprs <- dat.exprs
  
  dat.exprs.dup <- lapply(1:length(id.known), function(k){
    
    mod.dat.exprs[, id.missing] <-  mod.dat.exprs[, id.known[k]]
    
    mod.dat.exprs
    
  })
  
  
  lm.res <- lapply(1:length(dat.exprs.dup), function(i){
    
    dat.expr <- dat.exprs.dup[[i]]
    
    res <- lapply(1:length(study.id), function(k){
      
      y <- dat.aac[, study.id[k]]
      y <- as.numeric(scale(as.numeric(y), center = T, scale = T))
      x <- dat.expr[ ,study.id[k]]
      
      fit <- lm(y ~ as.numeric(x) + group.tissue)
      coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
      r.squared <- summary(fit)$r.squared 
      df <- fit$df.residual
      coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
      rownames(coef.fit) <- NULL
      colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
      coef.fit
      
    })
    
    res <- do.call(rbind, res)
    df.id <- paste("n", res$df + 10, sep=" = ")
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
     file="/results/MetaNoIndept/Pan/one.dup.bayes.RData")

#############################################################
#############################################################
## two duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)

####################################################################
## design of study
####################################################################

ccl.dat <- pan.ccl.ccle
ccl.dat$tissueid <- factor(ccl.dat$tissueid)
group.id <- unique(ccl.dat$tissueid)

group.tissue <- lapply(1:length(group.id), function(k){
  
  ifelse(ccl.dat$tissueid == group.id[k], 1, 0)
  
})

group.tissue <- do.call(cbind, group.tissue)
colnames( group.tissue) <- group.id

####################################################################
## Independent analysis
####################################################################

load("/results/Meta/Pan/meta.pan.RData")

meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                          meta.effect.size.res$drug == drugs[l] &
                                          meta.effect.size.res$method.meta == "DL", ]

meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                           I2 = meta.res.indept$I2,
                           te.val.r = meta.res.indept$TE.random,
                           ci.low = meta.res.indept$lower.random,
                           ci.up = meta.res.indept$upper.random)

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
      
      mod.dat.exprs[, id.missing] <-  mod.dat.exprs[, id.known[k]]
      
      mod.dat.exprs
      
    })
    
    lm.res <- lapply(1:length(dat.exprs.dup), function(i){
      
      dat.expr <- dat.exprs.dup[[i]]
      
      res <- lapply(1:length(study.id), function(k){
        
        y <- dat.aac[, study.id[k]]
        y <- as.numeric(scale(y, center = TRUE, scale = TRUE))
        x <- dat.expr[ ,study.id[k]]
        fit <- lm(y ~ as.numeric(x) + group.tissue)
        coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
        r.squared <- summary(fit)$r.squared 
        df <- fit$df.residual
        coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
        rownames(coef.fit) <- NULL
        colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
        coef.fit
        
      })
      
      res <- do.call(rbind, res)
      df.id <- paste("n", res$df + 10, sep=" = ")
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
     file="/results/MetaNoIndept/Pan/two.dup.bayes.RData")

#############################################################
#############################################################
## three duplication
#############################################################
#############################################################

dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
####################################################################
## design of study
####################################################################
ccl.dat <- pan.ccl.ccle
ccl.dat$tissueid <- factor(ccl.dat$tissueid)
group.id <- unique(ccl.dat$tissueid)

group.tissue <- lapply(1:length(group.id), function(k){
  
  ifelse(ccl.dat$tissueid == group.id[k], 1, 0)
  
})

group.tissue <- do.call(cbind, group.tissue)
colnames( group.tissue) <- group.id

####################################################################
## Independent analysis
####################################################################
load("/results/Meta/Pan/meta.pan.RData")

meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                          meta.effect.size.res$drug == drugs[l] &
                                          meta.effect.size.res$method.meta == "DL", ]

meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                           I2 = meta.res.indept$I2,
                           te.val.r = meta.res.indept$TE.random,
                           ci.low = meta.res.indept$lower.random,
                           ci.up = meta.res.indept$upper.random)

  
id <- lapply(1:(length(study.id)-(1+1)), function(k){
    
    id1 <- study.id[k]
    id2 <- study.id[k+1]
    
    id.dup <-  sapply((k+(1+1)):length(study.id), function(i){
      
      cbind(id1, id2, study.id[i])
      
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
      
      mod.dat.exprs[, id.missing] <-  mod.dat.exprs[, id.known[k]]
      
      mod.dat.exprs
      
    })
    
    lm.res <- lapply(1:length(dat.exprs.dup), function(i){
      
      dat.expr <- dat.exprs.dup[[i]]
      
      res <- lapply(1:length(study.id), function(k){
        
        y <- dat.aac[, study.id[k]]
        y <- as.numeric(scale(y, center = TRUE, scale = TRUE))
        x <- dat.expr[ ,study.id[k]]
        fit <- lm(y ~ as.numeric(x) + group.tissue)
        coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
        r.squared <- summary(fit)$r.squared 
        df <- fit$df.residual
        coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
        rownames(coef.fit) <- NULL
        colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
        coef.fit
        
      })
      
      res <- do.call(rbind, res)
      df.id <- paste("n", res$df + 10, sep=" = ")
      rownames(res) <- paste(study.id, df.id, sep=", ")
      id.miss <- paste(id.missing[1], id.missing[2],
                       id.missing[3], sep="-")
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
     file="/results/MetaNoIndept/Pan/three.dup.bayes.RData")


#############################################################
#############################################################
## four duplication
#############################################################
#############################################################


dat.meta.dup <- lapply(1:length(drugs), function(l){
  
  res.gene.drug <- lapply(1:length(id.gene), function(t){
    
    sub.ccle  <- z.exprs.ccle.data[z.exprs.ccle.data$gene_name == id.gene[t], ][-1]
    sub.ctrpv  <- z.exprs.ctrpv.data[z.exprs.ctrpv.data$gene_name == id.gene[t], ][-1]
    sub.gcsi  <- z.exprs.gcsi.data[z.exprs.gcsi.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v1  <- z.exprs.gdsc.v1.data[z.exprs.gdsc.v1.data$gene_name == id.gene[t], ][-1]
    sub.gdsc.v2  <- z.exprs.gdsc.v2.data[z.exprs.gdsc.v2.data$gene_name == id.gene[t], ][-1]
    
    sub.aac.ccle <- mod.aac.ccle[ ,colnames(mod.aac.ccle) == drugs[l] ]
    sub.aac.gcsi <- mod.aac.gcsi[ ,colnames(mod.aac.gcsi) == drugs[l] ]
    sub.aac.ctrpv <- mod.aac.ctrpv[ ,colnames(mod.aac.ctrpv) == drugs[l] ]
    sub.aac.gdsc.v1 <- mod.aac.gdsc.v1[ ,colnames(mod.aac.gdsc.v1) == drugs[l] ]
    sub.aac.gdsc.v2 <- mod.aac.gdsc.v2[ ,colnames(mod.aac.gdsc.v2) == drugs[l] ]
    
    expr.ccle <- scale(as.numeric(sub.ccle), center = T, scale = T)
    expr.ctrpv <- scale(as.numeric(sub.ctrpv), center = T, scale = T)
    expr.gcsi <- scale(as.numeric(sub.gcsi), center = T, scale = T)
    expr.gdsc.v1 <- scale(as.numeric(sub.gdsc.v1), center = T, scale = T)
    expr.gdsc.v2 <- scale(as.numeric(sub.gdsc.v2), center = T, scale = T)
    
    dat.exprs <- data.frame(CCLE = as.numeric(expr.ccle[,1]), 
                            CTRP = as.numeric(expr.ctrpv[,1]), 
                            gCSI = as.numeric(expr.gcsi[,1]), 
                            GDSC1 = as.numeric(expr.gdsc.v1[,1]),
                            GDSC2 = as.numeric(expr.gdsc.v2[,1]))
    
    dat.aac <- data.frame(CCLE = as.numeric( sub.aac.ccle), 
                          CTRP = as.numeric( sub.aac.ctrpv), 
                          gCSI = as.numeric( sub.aac.gcsi), 
                          GDSC1 = as.numeric( sub.aac.gdsc.v1),
                          GDSC2 = as.numeric( sub.aac.gdsc.v2))
    
    
    rownames(dat.exprs) <- colnames(sub.ccle)
    study.id <- colnames(dat.exprs)
    
####################################################################
## design of study
####################################################################
ccl.dat <- pan.ccl.ccle
ccl.dat$tissueid <- factor(ccl.dat$tissueid)
group.id <- unique(ccl.dat$tissueid)

group.tissue <- lapply(1:length(group.id), function(k){
  
  ifelse(ccl.dat$tissueid == group.id[k], 1, 0)
  
})

group.tissue <- do.call(cbind, group.tissue)
colnames( group.tissue) <- group.id

####################################################################
## Independent analysis
####################################################################
load("/results/Meta/Pan/meta.pan.RData")

meta.res.indept <- meta.effect.size.res[meta.effect.size.res$gene.name == id.gene[t] &
                                          meta.effect.size.res$drug == drugs[l] &
                                          meta.effect.size.res$method.meta == "DL", ]

meta.indept <-  data.frame(id.missing = "Indept", id.duplicate = "Indept",
                           I2 = meta.res.indept$I2,
                           te.val.r = meta.res.indept$TE.random,
                           ci.low = meta.res.indept$lower.random,
                           ci.up = meta.res.indept$upper.random)

res.all <- lapply(1:length(study.id), function(m){
  
  id.known <- study.id[m]
  id.missing <- study.id[study.id != id.known]
  mod.dat.exprs <- dat.exprs
  
  dat.exprs.dup <- lapply(1:length(id.known), function(k){
    
    mod.dat.exprs[, id.missing] <-  mod.dat.exprs[, id.known[k]]
    
    mod.dat.exprs
    
  })
  
  
  lm.res <- lapply(1:length(dat.exprs.dup), function(i){
    
    dat.expr <- dat.exprs.dup[[i]]
    
    res <- lapply(1:length(study.id), function(k){
      
      y <- dat.aac[, study.id[k]]
      y <- as.numeric(scale(y, center = TRUE, scale = TRUE))
      x <- dat.expr[ ,study.id[k]]
      fit <- lm(y ~ as.numeric(x) + group.tissue)
      coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(x)",]
      r.squared <- summary(fit)$r.squared 
      df <- fit$df.residual
      coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df)
      rownames(coef.fit) <- NULL
      colnames(coef.fit) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df")
      coef.fit
      
    })
    
    res <- do.call(rbind, res)
    df.id <- paste("n", res$df + 9, sep=" = ")
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
     file="/results/MetaNoIndept/Pan/four.dup.bayes.RData")



