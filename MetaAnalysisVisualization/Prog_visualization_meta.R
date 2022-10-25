
source("/code/MetaAnalysis/Prog_MetaAnalysis.run.R")
source("/code/MetaAnalysisVisualization/Prog_upset_meta.R")

dir.create("/results/Meta")
dir.create("/results/Meta/Breast")
dir.create("/results/Meta/Pan")
####################################################################
## UpSet plot to compare different meta-analysis methods: breast data
####################################################################

load("/results/Meta/Breast/meta.breast.RData")

id <- c("Erlotinib", "Lapatinib", "Paclitaxel")

for(i in 1:length(id)){
  
  res.adj.Fisher <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.fisher")]
  res.adj.Stouffer <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.stouffer")]
  res.adj.acat <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.acat")]
  res.adj.hmp <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.hmp")]
  
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
  sig.acat <- res.adj.acat[res.adj.acat$padj.acat < 0.05, "gene.name"]
  sig.hmp <- res.adj.hmp[res.adj.hmp$padj.hmp < 0.05, "gene.name"]
  
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
                    Stouffer = sig.stouffer,
                    CCT = sig.acat,
                    HMP = sig.hmp
                    )
  
  
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
  
  col.id <- c("maroon","blue2", "darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3", "salmon3", "springgreen4", "ivory4")
  sets = c("DL", "SJ", "HS", "HE", "EB", "PM", "Fisher", "Stouffer", "CCT","HMP")
  
  #col.id <- c("maroon","blue2", "darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3", "salmon3")
  #sets = c("DL", "SJ", "HS", "HE", "EB", "PM", "Fisher", "Stouffer")
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
  res.adj.acat <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.acat")]
  res.adj.hmp <- meta.pvalue.res[meta.pvalue.res$drug == id[i], c("gene.name", "padj.hmp")]
  
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
  sig.acat <- res.adj.acat[res.adj.acat$padj.acat < 0.05, "gene.name"]
  sig.hmp <- res.adj.hmp[res.adj.hmp$padj.hmp < 0.05, "gene.name"]
  
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
                    Stouffer = sig.stouffer,
                    CCT = sig.acat,
                    HMP = sig.hmp
                    )
  
  
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
  
  col.id <- c("maroon","blue2", "darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3", "salmon3", "springgreen4", "ivory4")
  sets = c("DL", "SJ", "HS", "HE", "EB", "PM", "Fisher", "Stouffer", "CCT","HMP")
  
  #col.id <- c("maroon","blue2", "darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3", "salmon3")
  #sets = c("DL", "SJ", "HS", "HE", "EB", "PM", "Fisher", "Stouffer")
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
  20, #length
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
    20, #length
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
    10, #length
    10, #width
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
    10, #length
    10, #width
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
  
  res.acat <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.acat <- res.acat[order(res.acat$acat.p), ]
  sub.acat <- res.acat[1:k, c("gene.name", "drug")]
  sub.acat <- data.frame(method = "ACAT", sub.acat)
  sub.acat$method <- "CCT"
  
  res.hmp <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.hmp <- res.hmp[order(res.hmp$hmp.p), ]
  sub.hmp <- res.hmp[1:k, c("gene.name", "drug")]
  sub.hmp <- data.frame(method = "HMP", sub.hmp)
  
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
  
  data <- rbind(sub.fisher, sub.stouffer, sub.acat, sub.hmp, sub.SJ, 
                sub.DL, sub.PM, sub.HS, sub.HE, sub.EB) #sub.acat, sub.hmp,
  
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
  
  res.acat <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.acat <- res.acat[order(res.acat$acat.p), ]
  sub.acat <- res.acat[1:k, c("gene.name", "drug")]
  sub.acat <- data.frame(method = "ACAT", sub.acat)
  sub.acat$method <- "CCT"
  
  res.hmp <- meta.pvalue.res[meta.pvalue.res$drug == id.drug[j], ]
  res.hmp <- res.hmp[order(res.hmp$hmp.p), ]
  sub.hmp <- res.hmp[1:k, c("gene.name", "drug")]
  sub.hmp <- data.frame(method = "HMP", sub.hmp)
  
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
  
  data <- rbind(sub.fisher, sub.stouffer, sub.acat, sub.hmp,sub.SJ, 
                sub.DL, sub.PM, sub.HS, sub.HE, sub.EB) # sub.acat, sub.hmp,  
  
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
    theme(axis.text.x=element_text(size=12,  face="bold"),
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
    theme(axis.text.x=element_text(size=12,  face="bold"),
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
  
  print(p)
  
  dev.off()
