
source("/code/MetaAnalysisVisualization/Prog_upset_meta.R")
source("/code/MetaAnalysis/Prog_MetaAnalysisAllCCLs.run.R")

dir.create("/results/Meta")
dir.create("/results/Meta/Breast/AllCCLs")

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
  
  jpeg(dir, res=600, width=3500, height=3000)
  
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



