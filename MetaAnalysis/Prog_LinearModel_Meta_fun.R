## Packages

library(meta)
#library(metap)
library(metafor)
library(GGally)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(RColorBrewer)
library(ggrepel)
library(bayesmeta)
library(forestplot)
library(trend)
################################################################################################
################################################################################################
############################### (adjusted) linear model functions ##############################
################################################################################################
################################################################################################

## Linear model fit to assess association between gene and drug sensitivity

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

## Adjusted linear model fit to assess association between gene and drug sensitivity

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

################################################################################################
################################################################################################
############################### meta-analysis functions ########################################
################################################################################################
################################################################################################

# Combining p-values: Fisher and Stouffer method

sumLog <- function (p){
  keep <- (p > 0) & (p <= 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(chisq = NA_real_, df = NA_integer_, p = NA_real_, 
                validp = p[keep])
  }
  else {
    lnp <- log(p[keep])
    chisq <- (-2) * sum(lnp)
    df <- 2 * length(lnp)
    if (length(lnp) != length(p)) {
      warning("Some studies omitted")
    }
    res <- list(chisq = chisq, df = df, p = pchisq(chisq, 
                                                   df, lower.tail = FALSE), validp = p[keep])
  }
  class(res) <- c("sumlog", "metap")
  res
}

sumZ <- function (p, weights = NULL, data = NULL, subset = NULL, na.action = na.fail){
  if (is.null(data)) 
    data <- sys.frame(sys.parent())
  mf <- match.call()
  mf$data <- NULL
  mf$subset <- NULL
  mf$na.action <- NULL
  mf[[1]] <- as.name("data.frame")
  mf <- eval(mf, data)
  if (!is.null(subset)) 
    mf <- mf[subset, ]
  mf <- na.action(mf)
  p <- as.numeric(mf$p)
  weights <- mf$weights
  noweights <- is.null(weights)
  if (noweights) 
    weights <- rep(1, length(p))
  if (length(p) != length(weights)) 
    warning("Length of p and weights differ")
  keep <- (p > 0) & (p < 1)
  invalid <- sum(1L * keep) < 2
  if (invalid) {
    warning("Must have at least two valid p values")
    res <- list(z = NA_real_, p = NA_real_, validp = p[keep], 
                weights = weights)
  }
  else {
    if (sum(1L * keep) != length(p)) {
      warning("Some studies omitted")
      omitw <- weights[!keep]
      if ((sum(1L * omitw) > 0) & !noweights) 
        warning("Weights omitted too")
    }
    zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep])/sqrt(sum(weights[keep]^2))
    res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE), 
                validp = p[keep], weights = weights)
  }
  class(res) <- c("sumz", "metap")
  res
}


fisher.fun <- function(pval){
  
  tmp <- sumLog(pval)
  tmp$p
  
}


stouffer.fun <- function(pval){
  
  tmp <- sumZ(pval)
  tmp$p
  
}


## combining effect sizes

meta.fun <- function(TE, seTE, method, study){
  
  res <- metagen(TE, seTE,  method.tau = method, studlab = study)
  
  res.val <- data.frame(TE.fixed = res$TE.fixed, seTE.fixed = res$seTE.fixed, 
                        lower.fixed = res$lower.fixed, upper.fixed = res$upper.fixed,
                        z.fixed = res$zval.fixed, pval.fixed = res$pval.fixed,
                        TE.random = res$TE.random, seTE.random = res$seTE.random, 
                        lower.random = res$lower.random, upper.random = res$upper.random,
                        z.random = res$zval.random, pval.random = res$pval.random,
                        I2 = res$I2, lower.I2 = res$lower.I2, upper.I2 = res$upper.I2, 
                        pval.Q = res$pval.Q, tau2 = res$tau2)
  res.val
  
}

## Bayesian meta-analysis

meta.bayes.fun <- function(TE, seTE, method, study){
  
  res <- bayesmeta(TE, 
                   seTE,
                   label=study,
                   mu.prior = c("mean"=NA, "sd"=NA), 
                   tau.prior=method)
  
  coef.res  <- res$summary["mean", "mu"]
  ci.res <- res$post.interval(mu.level=0.95, method="central")[c(1,2)]
  i2 <- round(res$I2(tau = res$summary["mean", "tau"]), 2) * 100
  
  res.val <- data.frame(coef = coef.res, ci.low = ci.res[1], 
                        ci.up = ci.res[2], I2 = i2)
  
  res.val
  
}


