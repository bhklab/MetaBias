## Packages

source("/code/Prog_LinearModel_Meta_fun.R")
source("/code/Prog_upset_meta.R")

dir.create("/results/LM")
dir.create("/results/LM/Breast")
dir.create("/results/LM/Breast/AllCCLs")

dir.create("/results/Meta")
dir.create("/results/Meta/Breast")
dir.create("/results/Meta/Breast/AllCCLs")

################################################################################################
################################################################################################
############################### (adjusted) linear model fitting ################################
################################################################################################
################################################################################################

## linear model: breast cancer cell line

load("/data/breast/DrugMeta.breast.z.RData")
load("/data/breast/DrugMeta.breast.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

z.aac <- list(mod.aac.ccle, mod.aac.ctrpv, mod.aac.gcsi,  mod.aac.gray, mod.aac.gdsc.v1, mod.aac.gdsc.v2, mod.aac.uhn)

z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data, z.exprs.gray.data, z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data, z.exprs.uhn.data)

study <- c("ccle", "ctrpv", "gcsi", "gray", "gdsc1", "gdsc2", "uhn")

for(i in 1:length(study)){
  
  print(i)
  lm.res <- lm.fun(exprs = z.exprs[[i]], aac= z.aac[[i]], drug = drugs)   
  save(lm.res, file = paste("/results/LM/Breast/AllCCLs", 
                                     paste("lm.breast", study[i], "RData", sep="."), sep="/") ) 
  
}

################################################################################################
################################################################################################
############################### meta-analysis fitting ##########################################
################################################################################################
################################################################################################

## Load Breast linear model results 

load("/results/LM/Breast/AllCCLs/lm.breast.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gray.RData")
res.gray <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.uhn.RData")
res.uhn <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.ctrpv.RData")
res.ctrpv <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gdsc2.RData")
res.gdsc2 <- lm.res


## merge data

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel") 
lm.dat <- list(res.ccle, res.ctrpv, res.gcsi, res.gray, res.gdsc1, res.gdsc2, res.uhn)
study <- c("ccle", "ctrp", "gcsi", "gray", "gdsc1", "gdsc2", "uhn")
method.meta <- c("DL")

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
     file = "/results/Meta/Breast/AllCCLs/meta.breast.RData")


#############################################################################################################
#############################################################################################################
##########################################  meta-analysis ###################################################
#############################################################################################################
#############################################################################################################

####################################################################
## forest plot: breast data
####################################################################
load("/results/LM/Breast/AllCCLs/lm.breast.ccle.RData")
res.ccle <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gcsi.RData")
res.gcsi <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gray.RData")
res.gray <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.uhn.RData")
res.uhn <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.ctrpv.RData")
res.ctrp <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gdsc1.RData")
res.gdsc1 <- lm.res
load("/results/LM/Breast/AllCCLs/lm.breast.gdsc2.RData")
res.gdsc2 <- lm.res

int <- c("Erlotinib")
id <- c("EGFR")

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
  
  
  dir0 <- "/results/Meta/Breast/AllCCLs"
  dir <- paste(paste(dir0, paste("SFig5_forest", 
                                 paste(int[k], id[k], sep="_"), sep="_"), sep="/"), "jpeg", sep=".")
  
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

load("/results/Meta/Breast/AllCCLs/meta.breast.RData")

id.drug <- c("Erlotinib")
id.gene <- c("EGFR")

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
  
  
  dir0 <- "/results/Meta/Breast/AllCCLs"
  dir <- paste(paste(dir0, paste("SFig5_volcano", paste(id.drug[i], id.gene[i], sep="_"), sep="_"), sep="/"), "jpeg", sep=".")
  
  jpeg(dir)
  
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
    geom_point(size = 1.2) + theme_minimal() +
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



