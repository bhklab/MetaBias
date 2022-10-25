
library(GGally)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(jaccard)
library(trend)
library(Cairo)

source("/code/MetaAnalysisNoIndept/Breast/Prog_duplication.R")
source("/code/MetaAnalysisNoIndept/Pan/Prog_duplication.R")
source("/code/MetaAnalysisNoIndept/Prog_correlation.run.R")

dir.create("/results/MetaDup")
dir.create("/results/MetaDup/Breast")
dir.create("/results/MetaDup/Pan")
####################################################################
## Breast : duplication
####################################################################
load("/results/MetaNoIndept/Breast/one.dup.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Breast/two.dup.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Breast/three.dup.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Breast/four.dup.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

load("/results/MetaNoIndept/Breast/five.dup.RData")
dat5 <- data.frame(missing.GEx = "five", dat.meta.dup)

load("/results/MetaNoIndept/Breast/six.dup.RData")
dat6 <- data.frame(missing.GEx = "six", dat.meta.dup)

load("/results/Cor/cor.breast.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id <- unique(dat1$gene_name)

res.gene.drug.all <- lapply(1:length(drugs), function(l){
  print(l)
  res.gene.drug <- lapply(1:length(id), function(k){
    print(k)
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
    
    sub.data$bias.te.f <- sub.data$te.val.f - indept.data$te.val.f[1]
    sub.data$bias.te.r <- sub.data$te.val.r - indept.data$te.val.r[1]
    
    data <- sub.data[sub.data$id.duplicate != "Indept", ]
    sub.data$id.duplicate <- factor(sub.data$id.duplicate)
    id.miss <- unique(data$missing.GEx)
    id.samples <- unique(data$samples)
    
    dat <- data[data$gene == id[k], ]
    
    mean.dat <- lapply(1:length(id.miss), function(j){
      
      bias.te.f.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.f))
      sd.te.f.val <- sd(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.f))
      
      bias.te.r.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))
      sd.te.r.val <- sd(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))  
      
      tau.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$tau)
      I2.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$I2)
      
      pval.Q.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$pval.Q)
      
      data.frame(id = id.miss[j], samples = id.samples[j], 
                 drug = unique(sub.data$drug), gene = id[k], 
                 I2.indept = unique(indept.data$I2),
                 pval.Q = unique(indept.data$pval.Q),
                 tau.mean = tau.val,
                 I2.mean = I2.val,
                 pval.Q.mean = pval.Q.val,
                 bias.te.f.val, bias.te.r.val,
                 sd.te.f.val, sd.te.r.val)
      
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
     file="/results/MetaNoIndept/Breast/cor.bias.res.RData")

#########################################################
##  Breast: Top genes to plot
#########################################################

load("/results/MetaNoIndept/Breast/cor.bias.res.RData")
load("/results/Meta/Breast/meta.breast.RData")

drugs <- "Erlotinib"
id <- c("CHMP7", "FBXO3", "DUOX2", "EGFR", "CD63", "ZNF143", 
        "PSMB3",  "C3orf52", "PITX3", "ERBB2", "CD84", "IGF2BP3", 
        "TGM3", "ABCG4", "S100A1", "COL11A1", "SCRT1")

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
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(2,5,3,4,1,6)])

p <-ggplot(data.plot, 
           aes(y = as.numeric(bias.te.r.val), x = id, fill= I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
  file = "/results/MetaDup/Breast/Fig5_Erlotinib.jpeg",
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
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(6,1,5,4,2,3)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  ylab("MAD pooled effect size, Lapatinib") +
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
  file = "/results/MetaDup/Breast/Fig5_Lapatinib.jpeg",
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
                                 res.gene.drug.all$drug == drugs , ], 
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
  ylab("MAD pooled effect size, Paclitaxel") +
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
  file = "/results/MetaDup/Breast/Fig5_Paclitaxel.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

###################################################
##  Breast: Correlation vs bias
##################################################

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
dat <- lapply(1:length(drugs), function(l){
  
  sub.dat <- cor.bias.res.all[cor.bias.res.all$drug == drugs[l], ]
  sub.dat$label <- NA
  sub.dat$label <- sapply(1:nrow(sub.dat), function(k){
    
    cor.val <- sub.dat$median.cor[k]
    if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
      if(abs(cor.val) > 0.7){  label.val <- "High"}else{
        label.val <- "Medium"
      }
    }
    
    label.val
    
  })
  
  sub.dat$label <- factor(sub.dat$label)
  sub.dat$label <- factor(sub.dat$label, levels = levels(sub.dat$label)[c(2,3,1)])
  sub.dat
  
  
})

dat <- do.call(rbind, dat)

p <-ggplot(dat, aes(label, bias.mean, fill = label)) +
  geom_violin(width = 0.8) +
  geom_boxplot(width=0.08) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  scale_fill_manual(values=c("indianred2", "cadetblue3", "#E69F00")) + 
  facet_wrap( ~ drug, ncol=3) +
  xlab("") +
  ylab("average of MAD of pooled effect size") +
  theme(axis.text.x=element_text(size=14,  face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=14),
        strip.text= element_text(size=12, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE)

Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaDup/Breast/Fig5_cor_bias.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)
print(p)

dev.off()

################################################################
## Breast: Jaccard Index
################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}

load("/results/Meta/Breast/meta.breast.RData")

load("/results/MetaNoIndept/Breast/one.dup.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Breast/two.dup.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Breast/three.dup.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Breast/four.dup.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

load("/results/MetaNoIndept/Breast/five.dup.RData")
dat5 <- data.frame(missing.GEx = "five", dat.meta.dup)

load("/results/MetaNoIndept/Breast/six.dup.RData")
dat6 <- data.frame(missing.GEx = "six", dat.meta.dup)

dat.file <- rbind(dat1, dat2, dat3, dat4, dat5, dat6) 
drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

res.JC <- lapply(1:length(drugs), function(l){
  
dat.indept <- meta.effect.size.res[meta.effect.size.res$drug == drugs[l]
                                   & meta.effect.size.res$method.meta == "DL", ]

dat.indept <- dat.indept[order(dat.indept$pval.random), ]
top.sig.gene <- dat.indept$gene.name[1:100]

data <- dat.file[dat.file$drug == drugs[l], ]
id.missing.GEx <- unique(data$missing.GEx)

sub.data <-  lapply(1:length(id.missing.GEx), function(i){

  dat <- data[data$missing.GEx == id.missing.GEx[i], ]
  
  id.missing <- unique(dat$id.missing)[-1]
  
  mod.dat <- lapply(1:length(id.missing), function(k){
    
    sub.dat <- dat[dat$id.missing == id.missing[k], ]
    id.duplicate <- unique(sub.dat$id.duplicate)
    
    dat.dup <- lapply(1:length(id.duplicate), function(j){
      
      sub.dat.dup <- sub.dat[sub.dat$id.duplicate == id.duplicate[j], ]
      sub.dat.dup[order(sub.dat.dup$pval.r), ][1:100, ]
      
    })
    
    do.call(rbind, dat.dup)
    
  })
  
  do.call(rbind, mod.dat)
  
})

sub.dat <- do.call(rbind, sub.data)


JI.similarity.res <- lapply(1:length(id.missing.GEx), function(i){
  
  id.missing <- unique(sub.dat[sub.dat$missing.GEx == id.missing.GEx[i], "id.missing"])
  
  similarity.res <- lapply(1:length(id.missing), function(k){
    
    id <- unique(sub.dat[sub.dat$id.missing == id.missing[k] &
                           sub.dat$missing.GEx == id.missing.GEx[i], "id.duplicate"])
    
    index.val <- sapply(1:length(id), function(j){
      
      jaccard(sub.dat[sub.dat$missing.GEx == id.missing.GEx[i] &
                        sub.dat$id.duplicate == id[j] & 
                        sub.dat$id.missing == id.missing[k], "gene_name"],
              dat.indept$gene.name[1:100])
      
    }) 
    
    data.frame(missing.GEx = id.missing.GEx[i],
               drug = sub.dat$drug[1],
               id.missing = id.missing[k],
               id.duplicate = id, 
               JI = index.val)
    
  })
  
  do.call(rbind,  similarity.res)
  
})

 do.call(rbind, JI.similarity.res)
 
 
 })

dat.plot <- do.call(rbind, res.JC)
dat.plot$missing.GEx <- factor(dat.plot$missing.GEx)
dat.plot$missing.GEx <- factor(dat.plot$missing.GEx,
                               levels = levels(dat.plot$missing.GEx)[c(3,6,5,2,1,4)])

p <-ggplot(dat.plot , aes(x= missing.GEx, JI, fill = missing.GEx)) +
  geom_violin(width = 0.8) +
  geom_boxplot(width=0.06) +
  stat_summary(fun.y=median, geom="point", size=1.5, color="black") +
  scale_fill_brewer(palette="Blues")  +
  facet_wrap( ~ drug, ncol=1) +
  xlab("") +
  ylab("Jaccard similarity index") +
  theme(axis.text.x=element_text(size=14,  face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12),
        strip.text= element_text(size=12, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE)


Cairo::Cairo(
  15, #length
  15, #width
  file = "/results/MetaDup/Breast/Fig7A.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

###############################################################
## Breast: Trend test
###############################################################

load("/results/MetaNoIndept/Breast/cor.bias.res.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

trend.res <- lapply(1:length(drugs), function(l){
  
  dat <- res.gene.drug.all[res.gene.drug.all$drug == drugs[l], ]
  id <- unique(dat$gene)
  
  trend.pval <- lapply(1:length(id), function(k){
    
    sub.dat <- dat[dat$gene == id[k], ]
    trend.fit.sign <- cs.test(sub.dat$bias.te.r.val)
    trend.fit.kendal <- mk.test(sub.dat$bias.te.r.val, alternative = "greater")
    data.frame(gene = id[k], drug = drugs[l], 
               pval.sign = trend.fit.sign$p.value, pval.kendal = trend.fit.kendal$p.value)
    
  })
  
  do.call(rbind, trend.pval)
  
})

trend.res <- do.call(rbind,  trend.res)

write.csv(trend.res, file="/results/MetaDup/Breast/trend.test.csv")


#######################################################################
#######################################################################
## Pan data
####################################################################### 
#######################################################################

####################################################################
## Pan : duplication
####################################################################
load("/results/MetaNoIndept/Pan/one.dup.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Pan/two.dup.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Pan/three.dup.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Pan/four.dup.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

load("/results/Cor/cor.pan.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
id <- unique(dat1$gene_name)

res.gene.drug.all <- lapply(1:length(drugs), function(l){
  
  print(l)
  res.gene.drug <- lapply(1:length(id), function(k){
    
    print(k)
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
    
    sub.data$bias.te.f <- sub.data$te.val.f - indept.data$te.val.f[1]
    sub.data$bias.te.r <- sub.data$te.val.r - indept.data$te.val.r[1]
    
    data <- sub.data[sub.data$id.duplicate != "Indept", ]
    sub.data$id.duplicate <- factor(sub.data$id.duplicate)
    id.miss <- unique(data$missing.GEx)
    id.samples <- unique(data$samples)
    
    dat <- data[data$gene == id[k], ]
    
    mean.dat <- lapply(1:length(id.miss), function(j){
      
      bias.te.f.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.f))
      sd.te.f.val <- sd(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.f))
      
      bias.te.r.val <- mean(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))
      sd.te.r.val <- sd(abs(dat[dat$missing.GEx == id.miss[j], ]$bias.te.r))  
      
      tau.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$tau)
      I2.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$I2)
      
      pval.Q.val <- mean(dat[dat$missing.GEx == id.miss[j], ]$pval.Q)
      
      data.frame(id = id.miss[j], samples = id.samples[j], 
                 drug = unique(sub.data$drug), gene = id[k], 
                 I2.indept = unique(indept.data$I2),
                 pval.Q = unique(indept.data$pval.Q),
                 tau.mean = tau.val,
                 I2.mean = I2.val,
                 pval.Q.mean = pval.Q.val,
                 bias.te.f.val, bias.te.r.val,
                 sd.te.f.val, sd.te.r.val)
      
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
     file="/results/MetaNoIndept/Pan/cor.bias.res.RData")

#########################################################
##  Pan: Top genes to plot
#########################################################

load("/results/MetaNoIndept/Pan/cor.bias.res.RData")
load("/results/Meta/Pan/meta.pan.RData")

drugs <- "Erlotinib"
id <- c("C1orf116", "SPRR1B", "KRT1", "EGFR", "COX7A2", "ROS1" ,
        "TDRD1", "ZNF365", "KRT1", "ERBB2", "MAPK7", "CYP2A13", 
        "C4orf19", "ZNF688", "CPSF6", "S100A1", "NOC3L", "CD244")

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
dat.low.i2 <- data.frame(dat.low.i2, I2.id = "I2 < 50%")

dat.plot <- rbind(dat.high.i2, dat.low.i2)
data.plot <- lapply(1:length(dat.plot$gene.name), function(k){
  
  data.frame(res.gene.drug.all[res.gene.drug.all$gene ==  dat.plot$gene.name[k] &
                                 res.gene.drug.all $drug == drugs , ], 
             dat.plot[dat.plot$gene.name == dat.plot$gene.name[k], -c(1,2)])
})

data.plot <- do.call(rbind, data.plot)
data.plot$id <- factor(data.plot$id)
data.plot$id <- factor(data.plot$id, levels = levels(data.plot$id)[c(2,4,3,1)])
data.plot$gene <- factor(data.plot$gene)

cor.id <- paste("r = ", round(data.plot$median.cor, 2), sep="")
data.plot$o.gene <- paste(data.plot$gene, cor.id, sep=", ")
data.plot$o.gene <- factor(data.plot$o.gene)
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(1, 6, 4, 3,2,5)])

p <-ggplot(data.plot, 
           aes(y = as.numeric(bias.te.r.val), x = id, fill= I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
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
  file = "/results/MetaDup/Pan/Fig6_Erlotinib.jpeg",
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
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(5, 6, 3, 2, 4, 1)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  ylab("MAD pooled effect size, Lapatinib") +
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
  file = "/results/MetaDup/Pan/Fig6_Lapatinib.jpeg",
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
data.plot$o.gene <- factor(data.plot$o.gene, levels = levels(data.plot$o.gene)[c(1, 6,3,5,4,2)])

p <-ggplot(data.plot , aes(y = bias.te.r.val, x = id, fill=I2.id) )+
  geom_bar(width = 0.4, stat="identity", colour = "grey20") +
  facet_wrap(. ~ o.gene) +
  scale_fill_manual(values=c("#999999", "#56B4E9")) +
  ylab("MAD pooled effect size, Paclitaxel") +
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
  file = "/results/MetaDup/Pan/Fig6_Paclitaxel.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()


###################################################
##  Pan: Correlation vs bias
##################################################

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")
dat <- lapply(1:length(drugs), function(l){
  
  sub.dat <- cor.bias.res.all[cor.bias.res.all$drug == drugs[l], ]
  sub.dat$label <- NA
  sub.dat$label <- sapply(1:nrow(sub.dat), function(k){
    
    cor.val <- sub.dat$median.cor[k]
    if(abs(cor.val) < 0.3){ label.val <- "Low"  }else{
      if(abs(cor.val) > 0.7){  label.val <- "High"}else{
        label.val <- "Medium"
      }
    }
    
    label.val
    
  })
  
  sub.dat$label <- factor(sub.dat$label)
  sub.dat$label <- factor(sub.dat$label, levels = levels(sub.dat$label)[c(2,3,1)])
  sub.dat
  
  
})

dat <- do.call(rbind, dat)

p <-ggplot(dat, aes(label, bias.mean, fill = label)) +
  geom_violin(width = 0.8) +
  geom_boxplot(width=0.08) +
  stat_summary(fun.y=median, geom="point", size=2, color="black") +
  scale_fill_manual(values=c("indianred2", "cadetblue3", "#E69F00")) + 
  facet_wrap( ~ drug, ncol=3) +
  xlab("") +
  ylab("average of MAD of pooled effect size") +
  theme(axis.text.x=element_text(size=14,  face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=14),
        strip.text= element_text(size=12, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE)


Cairo::Cairo(
  20, #length
  15, #width
  file = "/results/MetaDup/Pan/Fig6_cor_bias.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

################################################################
## Pan: Jaccard Index
################################################################
jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}


load("/results/Meta/Pan/meta.pan.RData")

load("/results/MetaNoIndept/Pan/one.dup.RData")
dat1 <- data.frame(missing.GEx = "one", dat.meta.dup)

load("/results/MetaNoIndept/Pan/two.dup.RData")
dat2 <- data.frame(missing.GEx = "two", dat.meta.dup)

load("/results/MetaNoIndept/Pan/three.dup.RData")
dat3 <- data.frame(missing.GEx = "three", dat.meta.dup)

load("/results/MetaNoIndept/Pan/four.dup.RData")
dat4 <- data.frame(missing.GEx = "four", dat.meta.dup)

dat.file <- rbind(dat1, dat2, dat3, dat4) 
drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

res.JC <- lapply(1:length(drugs), function(l){
  
  dat.indept <- meta.effect.size.res[meta.effect.size.res$drug == drugs[l]
                                     & meta.effect.size.res$method.meta == "DL", ]
  
  dat.indept <- dat.indept[order(dat.indept$pval.random), ]
  top.sig.gene <- dat.indept$gene.name[1:100]
  
  data <- dat.file[dat.file$drug == drugs[l], ]
  id.missing.GEx <- unique(data$missing.GEx)
  
  sub.data <-  lapply(1:length(id.missing.GEx), function(i){
    
    dat <- data[data$missing.GEx == id.missing.GEx[i], ]
    
    id.missing <- unique(dat$id.missing)[-1]
    
    mod.dat <- lapply(1:length(id.missing), function(k){
      
      sub.dat <- dat[dat$id.missing == id.missing[k], ]
      id.duplicate <- unique(sub.dat$id.duplicate)
      
      dat.dup <- lapply(1:length(id.duplicate), function(j){
        
        sub.dat.dup <- sub.dat[sub.dat$id.duplicate == id.duplicate[j], ]
        sub.dat.dup[order(sub.dat.dup$pval.r), ][1:100, ]
        
      })
      
      do.call(rbind, dat.dup)
      
    })
    
    do.call(rbind, mod.dat)
    
  })
  
  sub.dat <- do.call(rbind, sub.data)
  
  
  JI.similarity.res <- lapply(1:length(id.missing.GEx), function(i){
    
    id.missing <- unique(sub.dat[sub.dat$missing.GEx == id.missing.GEx[i], "id.missing"])
    
    similarity.res <- lapply(1:length(id.missing), function(k){
      
      id <- unique(sub.dat[sub.dat$id.missing == id.missing[k] &
                             sub.dat$missing.GEx == id.missing.GEx[i], "id.duplicate"])
      
      index.val <- sapply(1:length(id), function(j){
        
        jaccard(sub.dat[sub.dat$missing.GEx == id.missing.GEx[i] &
                          sub.dat$id.duplicate == id[j] & 
                          sub.dat$id.missing == id.missing[k], "gene_name"],
                dat.indept$gene.name[1:100])
        
      }) 
      
      data.frame(missing.GEx = id.missing.GEx[i],
                 drug = sub.dat$drug[1],
                 id.missing = id.missing[k],
                 id.duplicate = id, 
                 JI = index.val)
      
    })
    
    do.call(rbind,  similarity.res)
    
  })
  
  do.call(rbind, JI.similarity.res)
  
  
})

dat.plot <- do.call(rbind, res.JC)
dat.plot$missing.GEx <- factor(dat.plot$missing.GEx)
dat.plot$missing.GEx <- factor(dat.plot$missing.GEx,
                               levels = levels(dat.plot$missing.GEx)[c(2,4,3,1)])

p <-ggplot(dat.plot , aes(x= missing.GEx, JI, fill = missing.GEx)) +
  geom_violin(width = 0.8) +
  geom_boxplot(width=0.06) +
  stat_summary(fun.y=median, geom="point", size=1.5, color="black") +
  scale_fill_brewer(palette="Blues")  +
  facet_wrap( ~ drug, ncol=1) +
  xlab("") +
  ylab("Jaccard similarity index") +
  theme(axis.text.x=element_text(size=14,  face="bold"),
        axis.title=element_text(size=14, face="bold"),
        axis.text.y=element_text(size=12),
        strip.text= element_text(size=12, face="bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  guides(fill=FALSE)

Cairo::Cairo(
  15, #length
  15, #width
  file = "/results/MetaDup/Pan/Fig7B.jpeg",
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

print(p)

dev.off()

###############################################################
## Pan: Trend test
###############################################################

load("/results/MetaNoIndept/Pan/cor.bias.res.RData")

drugs <- c("Erlotinib", "Lapatinib", "Paclitaxel")

trend.res <- lapply(1:length(drugs), function(l){
  
  dat <- res.gene.drug.all[res.gene.drug.all$drug == drugs[l], ]
  id <- unique(dat$gene)
  
  trend.pval <- lapply(1:length(id), function(k){
    
    sub.dat <- dat[dat$gene == id[k], ]
    trend.fit.sign <- cs.test(sub.dat$bias.te.r.val)
    trend.fit.kendal <- mk.test(sub.dat$bias.te.r.val, alternative = "greater")
    data.frame(gene = id[k], drug = drugs[l], 
               pval.sign = trend.fit.sign$p.value, pval.kendal = trend.fit.kendal$p.value)
    
  })
  
  do.call(rbind, trend.pval)
  
})

trend.res <- do.call(rbind,  trend.res)

write.csv(trend.res, file="/results/MetaDup/Pan/trend.test.csv")

