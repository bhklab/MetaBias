## Packages

library(meta)
#library(metap)
library(metafor)
library(harmonicmeanp)
library(ACAT)
library(GGally)
library(ggplot2)
library(dplyr)
library(UpSetR)
library(RColorBrewer)
library(ggrepel)
library(bayesmeta)
library(forestplot)
###################################################
## functions
###################################################

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

