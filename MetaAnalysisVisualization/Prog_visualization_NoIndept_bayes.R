
library(GGally)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)


source("/code/MetaAnalysisNoIndept/Breast/Prog_duplication_bayes.R")
source("/code/MetaAnalysisNoIndept/Pan/Prog_duplication_bayes.R")
source("/code/MetaAnalysisNoIndept/Prog_correlation.run.R")

dir.create("/results/MetaDup")
dir.create("/results/MetaDup/Breast")
dir.create("/results/MetaDup/Pan")
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
load("/results/Meta/Breast/meta.breast.RData")

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
  file = "/results/MetaDup/Breast/SFig12_Erlotinib.jpeg",
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
  file = "/results/MetaDup/Breast/SFig12_Lapatinib.jpeg",
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
  file = "/results/MetaDup/Breast/SFig12_Paclitaxel.jpeg",
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
  file = "/results/MetaDup/Breast/SFig12_CI_DL_Bayes.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

#######################################################################
#######################################################################
## Pan data
####################################################################### 
#######################################################################

####################################################################
## Pan : duplication
####################################################################
load("/results/MetaNoIndept/Pan/one.dup.bayes.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Pan/two.dup.bayes.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Pan/three.dup.bayes.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Pan/four.dup.bayes.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

load("/results/Cor/cor.pan.RData")

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
    
    sub.data <- rbind(sub.dat1, sub.dat2, sub.dat3,
                      sub.dat4)
    
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
  print(l)
  cor.bias.res <- lapply(1:length(id), function(k){
    print(k)
    sub.cor.res <- cor.res[cor.res$gene == id[k], ]
    sub.res <- res.gene.drug.all[res.gene.drug.all$gene == id[k] & 
                                   res.gene.drug.all$drug == drugs[l], ]
    I2.mean <- mean(sub.res$I2.mean)
    bias.mean <- mean(sub.res$bias.te.r.val)
    ci.mean <- mean(sub.res$ci.val)
    data.frame(gene = unique(sub.res$gene), 
               drug = unique(sub.res$drug), 
               bias.mean = bias.mean, 
               I2.mean = I2.mean, 
               I2 = sub.res$I2.indept[1],
               ci.mean = ci.mean, 
               sub.cor.res[, 2:3])
    
  })
  
  do.call(rbind, cor.bias.res)
  
})

cor.bias.res.all <- do.call(rbind, cor.bias.res.all)

save(cor.bias.res.all,  res.gene.drug.all ,
     file="/results/MetaNoIndept/Pan/cor.bias.bayes.res.RData")

#########################################################
##  Pan: Top genes to plot
#########################################################

load("/results/MetaNoIndept/Pan/cor.bias.bayes.res.RData")
load("/results/Meta/Pan/meta.pan.RData")

drugs <- "Erlotinib"
id <- c("C1orf116", "SPRR1B", "KRT1", "EGFR", "COX7A2", "ROS1" ,
        "TDRD1", "ZNF365", "KRT1", "ERBB2", "MAPK7", "CYP2A13", 
        "C4orf19", "ZNF688", "CPSF6", "S100A1", "NOC3L", "CD244")

gene.id <- id[1:6]
res <- meta.effect.size.res[meta.effect.size.res$drug == drugs & 
                              meta.effect.size.res$method.meta == "DL", ]
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
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all$drug ==drugs  , ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(2,4,3,1)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(1, 6, 4, 3, 2, 5)])

p <-ggplot(data.plot, 
           aes(y = as.numeric(bias.te.r.val), x = id, fill= I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  # ggtitle(plot.id) +
  ylab("MAD pooled effect size, Erlotinib") +
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
        legend.text = element_text(size = 14, face="bold"),
        legend.title = element_blank())

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaDup/Pan/SFig13_Erlotinib.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()


drugs <- "Lapatinib"
gene.id <- id[7:12]
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
                                 res.gene.drug.all$drug == drugs, ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(2,4,3,1)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(5, 6, 3, 2, 4, 1)])

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
  file = "/results/MetaDup/Pan/SFig13_Lapatinib.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

drugs <- "Paclitaxel"
gene.id <- id[13:18]
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
                                 res.gene.drug.all$drug == drugs, ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(2,4,3,1)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(1, 6,3,5,4,2)])

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
  file = "/results/MetaDup/Pan/SFig13_Paclitaxel.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)


print(p)

dev.off()

###################################################
##  Pan: 95% CI DL vs Bayes credible interval
##################################################

load("/results/MetaNoIndept/Pan/cor.bias.bayes.res.RData")

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
  
  load("/results/MetaNoIndept/Pan/one.dup.bayes.RData")
  dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)
  dat1 <- dat1[dat1$id.missing != "Indept", ]
  dat1 <- dat1[dat1$drug == drugs[l], ]
  dl.one <- mean(dat1$ci.up - dat1$ci.low)
  
  load("/results/MetaNoIndept/Pan/two.dup.bayes.RData")
  dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)
  dat2 <- dat2[dat2$id.missing != "Indept", ]
  dat2 <- dat2[dat2$drug == drugs[l], ]
  dl.two <- mean(dat2$ci.up - dat2$ci.low)
  
  load("/results/MetaNoIndept/Pan/three.dup.bayes.RData")
  dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)
  dat3 <- dat3[dat3$id.missing != "Indept", ]
  dat3 <- dat3[dat3$drug == drugs[l], ]
  dl.three <- mean(dat3$ci.up - dat3$ci.low)
  
  load("/results/MetaNoIndept/Pan/four.dup.bayes.RData")
  dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)
  dat4 <- dat4[dat4$id.missing != "Indept", ]
  dat4 <- dat4[dat4$drug == drugs[l], ]
  dl.four <- mean(dat4$ci.up - dat4$ci.low)
  
  data.frame(y = c(dl.one, dl.two, dl.three, dl.four,
                   bayes.one, bayes.two, bayes.three, bayes.four),
             drug = drugs[l],
             group = rep(c("one", "two", "three", "four"), 2), 
             method = c(rep("DerSimonian-Laird (DL)", 4), rep("Jeffreys", 4)))
})

dat <- do.call(rbind, data.plot)
dat$group <- factor(dat$group)
dat$group <- factor(dat$group, levels = levels(dat$group)[c(2,4,3,1)])

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
  file = "/results/MetaDup/Pan/SFig13_CI_DL_Bayes.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

