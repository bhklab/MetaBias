####################################################################
## Linear model fit to assess association between gene and drug sensitivity
####################################################################

lm.fun <- function(exprs, aac, drug){
  
  mod.aac <- t(aac)
  exprs.dat <- exprs
  
  res.rna.drug <- lapply(1:length(drug), function(k){
    
    mod.val <- mod.aac[rownames(mod.aac) == drug[k], ]
    mod <- (mod.val - mean(mod.val))/sd(mod.val)
    
    lm.fit <- lapply(1:nrow(exprs.dat), function(j){
      
      mod.dat <- exprs.dat[j,-1]
      fit <- lm(mod ~ as.numeric(mod.dat))
      coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(mod.dat)",]
      r.squared <- summary(fit)$r.squared 
      df <- fit$df.residual
      ci.res <- confint(fit)["as.numeric(mod.dat)", ]
      coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df, 
                             ci.low = ci.res[1], ci.high = ci.res[2])
      rownames(coef.fit) <- NULL
      coef.fit
      
    })
    
    coef.lm <- do.call(rbind, lm.fit)
    colnames(coef.lm) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df",
                           "ci.low", "ci.high")
    mod.coef.lm <- data.frame(ensembl.id = rownames(exprs.dat),
                              gene.name = exprs.dat$gene_name, 
                              drug = drug[k], coef.lm)
    mod.coef.lm
    
  })
  
  res.drug.rna <- do.call(rbind, res.rna.drug)
  res.rnaseq.drug <- res.drug.rna[!is.na(res.drug.rna$Pvalue), ]
  res.rnaseq.drug
  
}

####################################################################
## Adjusted linear model fit to assess association between gene and drug sensitivity
####################################################################

lm.adjusted.fun <- function(exprs, aac, drug, ccl){
  
  mod.aac <- t(aac)
  mod.aac <- mod.aac[, rownames(ccl)]
  exprs.dat <- exprs[, c("gene_name", colnames(mod.aac))]
  ccl$tissueid <- factor(ccl$tissueid)
  
  group.id <- unique(ccl$tissueid)
  
  group.tissue <- lapply(1:length(group.id), function(k){
    
    ifelse(ccl$tissueid == group.id[k], 1, 0)
    
  })
  
  group.tissue <- do.call(cbind, group.tissue)
  colnames( group.tissue) <- group.id
  
res.rna.drug <- lapply(1:length(drug), function(k){
  
  mod <- mod.aac[rownames(mod.aac) == drug[k], ]
  mod <- as.numeric(scale(mod, scale=TRUE, center = TRUE))
  
  
  lm.fit <- lapply(1:nrow(exprs.dat), function(j){
    
    mod.dat <- exprs.dat[j,-1]
    fit <- lm(mod ~ as.numeric(mod.dat) + group.tissue)
    coef.fit <- as.data.frame(summary(fit)$coefficients)["as.numeric(mod.dat)",]
    r.squared <- summary(fit)$r.squared 
    df <- fit$df.residual
    ci.res <- confint(fit)["as.numeric(mod.dat)", ]
    coef.fit <- data.frame(coef.fit, r2 = r.squared, df = df, 
                           ci.low = ci.res[1], ci.high = ci.res[2])
    rownames(coef.fit) <- NULL
    coef.fit
    
  })
  
  coef.lm <- do.call(rbind, lm.fit)
  colnames(coef.lm) <- c("estimate", "std", "tvalue", "Pvalue", "r2", "df",
                         "ci.low", "ci.high")
  mod.coef.lm <- data.frame(ensembl.id = rownames(exprs.dat),
                            gene.name = exprs.dat$gene_name, 
                            drug = drug[k], coef.lm)
  mod.coef.lm
  
})

res.drug.rna <- do.call(rbind, res.rna.drug)
res.rnaseq.drug <- res.drug.rna[!is.na(res.drug.rna$Pvalue), ]
res.rnaseq.drug 

}

