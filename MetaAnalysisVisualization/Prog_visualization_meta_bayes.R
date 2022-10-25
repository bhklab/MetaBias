source("/code/MetaAnalysis/Prog_MetaAnalysisBayes.run.R")
source("/code/MetaAnalysis/Prog_MetaAnalysis.run.R")

dir.create("/results/MetaVisualization")
dir.create("/results/MetaVisualization/Breast")
dir.create("/results/MetaVisualization/Pan")
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
    theme(axis.text.x=element_text(size=12,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
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
    theme(axis.text.x=element_text(size=12,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()
  
}


####################################################################
## Pan data: compare Bayesian and DL meta-analyses (I2 and overall effect size) 
####################################################################

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
    theme(axis.text.x=element_text(size=12,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
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
    theme(axis.text.x=element_text(size=12,  face="bold"),
          axis.title=element_text(size=10,face="bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black")) +
    guides(fill=FALSE)
  
  print(p)
  
  dev.off()

}

