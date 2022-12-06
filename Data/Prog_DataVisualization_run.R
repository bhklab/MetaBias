
source("/code/Prog_upset.R")

dir.create("/results/Summary_Figure")
dir.create("/results/Summary_Figure/Data")
#############################################################
## Fig1A: Breast Cancer Cell line
#############################################################
load("/data/breast/DrugMeta.breast.RData")

id.ccle <- rownames(brca.ccl.ccle)
id.ctrpv <- rownames(brca.ccl.ctrpv)
id.gcsi <- rownames(brca.ccl.gcsi)
id.gray <- rownames(brca.ccl.gray)
id.gdsc1 <- rownames(brca.ccl.gdsc.v1)
id.gdsc2 <- rownames(brca.ccl.gdsc.v2)
id.uhn <- rownames(brca.ccl.uhn)

listInput <- list(CCLE = id.ccle, 
                  CTRP = id.ctrpv,
                  gCSI = id.gcsi,
                  GRAY = id.gray,
                  GDSC1 = id.gdsc1, 
                  GDSC2 = id.gdsc2, 
                  UHNBreast = id.uhn)

col.id <- c("maroon","blue2","darkgoldenrod1", "darkorchid3", "deeppink1", "green2", "lightskyblue3")
sets = c("CCLE", "CTRP", "gCSI", "GRAY", "GDSC1", "GDSC2","UHNBreast")
mainbar.y.label = "breast cancer cell line";
sets.bar.color= col.id;

Cairo::Cairo(
  20, #length
  15, #width
  file = paste("/results/Summary_Figure/Data/Fig1A", ".jpeg", sep = ""),
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)

 upset.fun(fromList(listInput)) 

dev.off()

#############################################################
## Fig1B: upset plot Pan Cancer Cell line
#############################################################
load("/data/pan/DrugMeta.pan.mod.RData")

id.ccle <- rownames(pan.ccl.ccle)
id.ctrpv <- rownames(pan.ccl.ctrpv)
id.gcsi <- rownames(pan.ccl.gcsi)
id.gdsc1 <- rownames(pan.ccl.gdsc.v1)
id.gdsc2 <- rownames(pan.ccl.gdsc.v2)

listInput <- list(CCLE = id.ccle, 
                  CTRP = id.ctrpv,
                  gCSI = id.gcsi,
                  GDSC1 = id.gdsc1, 
                  GDSC2 = id.gdsc2)


col.id <- c("maroon","blue2","darkgoldenrod1", "deeppink1", "green2")
sets = c("CCLE", "CTRP", "gCSI", "GDSC1", "GDSC2")
sets.bar.color= col.id;
mainbar.y.label = "pan-cancer cell line";

 Cairo::Cairo(
   20, #length
   15, #width
   file = paste("/results/Summary_Figure/Data/Fig1B", ".jpeg", sep = ""),
   type = "jpeg", 
   bg = "transparent", 
   dpi = 300,
   units = "cm"
 )
 
upset.fun(fromList(listInput)) 

dev.off() 

###########################################
## Fig1B: pie chart Pan Cancer Cell line
###########################################

load("/data/pan/DrugMeta.pan.fixsample.mod.RData")

data <- data.frame(
  group=c("Bowel", "Breast", "Esophagus/Stomach", "Lung",
          "Ovary/Fallopian Tube", "Pancreas", "Skin"),
  value = as.numeric(table(pan.ccl.ccle$tissueid))
)

data <- data %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(data$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

Cairo::Cairo(
   10, #length
   15, #width
   file = paste("/results/Summary_Figure/Data/Fig1B_pie", ".jpeg", sep = ""),
   type = "jpeg", 
   bg = "transparent", 
   dpi = 300,
   units = "cm"
 )
 
ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="grey60") +
  coord_polar("y", start=0) +
  theme_void() + 
  geom_text(aes(y = ypos, label = paste(round(prop), "%", sep="")),
            color = "grey10", size=5) +
  theme(legend.title=element_blank(),
        legend.text = element_text(size = 6, face = "bold"),
        legend.position="top") +
  scale_fill_manual(values=c("indianred2", "dodgerblue3",
                             "darkorchid3", "chartreuse3",
                             "yellow3", "deeppink1",
                             "lightblue4")) 

dev.off()

##############################################################
## SFig2: bar plot Pan cancer data
#############################################################

load("/data/pan/DrugMeta.pan.mod.RData")

id.study <- c("GDSC1", "GDSC2", "CCLE", "CTRP", "gCSI")

pan.dat <- list(pan.ccl.gdsc.v1, pan.ccl.gdsc.v2, pan.ccl.ccle,
                pan.ccl.ctrpv, pan.ccl.gcsi)

id.dat <- names(table(pan.dat[[1]]$tissueid))
n.dat <- table(pan.dat[[1]]$tissueid)
names(n.dat) <- NULL
freq.dat <- round(n.dat/nrow(pan.dat[[1]]) * 100)
res.dat <- data.frame( id = id.dat, study = id.study[1], n.ccls = n.dat)[, -c(3)]

set.seed(23561)
n <- nrow(res.dat)
cpalette <- c("#D1E354", "#75D3E6", "#D7C2DD", "#8970D7", "#CDE5DA", "#799CD8", "#9B6F88", "#E45D8A", "#78E5C8", "#779D84", "#72E75B", "#81E797", "#E1B657", "#DBB39B", "#953CE5", "#DB56CC", "#D5E4A2", "#DB774E", "#DF98D9")
dat.col = data.frame(col = cpalette[1:length(id.dat)], id = id.dat)

dat <- rbind(res.dat)
dat$study <- factor(dat$study)
dat$x <- factor(dat$id,                                   
                levels = dat$id[order(dat$n.ccls.Freq, decreasing = FALSE)])
col.dat <- data.frame( id = cpalette, tissue = id.dat)
id.not <- NA
title.id <- paste("cancer cell lines", id.study[1], sep=", ")

p <-ggplot(dat , aes(x = x, y = n.ccls.Freq, fill = id)) +
  geom_bar(stat="identity", colour="black", width = 0.8) +
  coord_flip() +
  xlab(" ") +
  ylab(title.id) +
  theme(axis.text.x=element_text(size=10,  face="bold"),
        axis.title=element_text(size=12,face="bold"),
        axis.text.y=element_text(size=12, face = "bold"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position="none",
        legend.text = element_text(size = 9, face="bold"),
        legend.title = element_blank()) +
  scale_fill_manual(values=col.dat[!(col.dat$tissue %in% id.not), ]$id)

Cairo::Cairo(
  15, #length
  10, #width
  file = paste("/results/Summary_Figure/Data/SFig2_GDSC1", ".jpeg", sep = ""),
  type = "jpeg", 
  bg = "transparent", 
  dpi = 300,
  units = "cm"
)
  
p

dev.off()


id.dat0 <- names(table(pan.dat[[1]]$tissueid))

for(i in 2:length(id.study)){
  
  id.dat <- names(table(pan.dat[[i]]$tissueid))
  n.dat <- table(pan.dat[[i]]$tissueid)
  names(n.dat) <- NULL
  freq.dat <- round(n.dat/nrow(pan.dat[[i]]) * 100)
  res.dat <- data.frame( id = id.dat, study = id.study[i], n.ccls = n.dat)[, -c(3)]
  int <- intersect(id.dat, id.dat0)
  id.not <- setdiff(id.dat0, int)
  
  dat <- rbind(res.dat)
  dat$study <- factor(dat$study)
  dat$x <- factor(dat$id,                                   
                  levels = dat$id[order(dat$n.ccls.Freq, decreasing = FALSE)])
  col.dat <- data.frame(id = cpalette, tissue = id.dat0)
  
  title.id <- paste("cancer cell lines", id.study[i], sep=", ")
  
  p <-ggplot(dat , aes(x = x, y = n.ccls.Freq, fill = id)) +
    geom_bar(stat="identity", colour="black", width = 0.8) +
    coord_flip() +
    xlab(" ") +
    ylab(title.id) +
    theme(axis.text.x=element_text(size=10,  face="bold"),
          axis.title=element_text(size=12,face="bold"),
          axis.text.y=element_text(size=12, face = "bold"),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          plot.background = element_blank(), 
          axis.line = element_line(colour = "black"),
          legend.position="none",
          legend.text = element_text(size = 9, face="bold"),
          legend.title = element_blank()) +
    scale_fill_manual(values=col.dat[!(col.dat$tissue %in% id.not), ]$id)

  dir <- paste("/results/Summary_Figure/Data", paste(paste("SFig2", id.study[i], sep="_"), "jpeg", sep="."), sep="/")
  
  Cairo::Cairo(
    15, #length
    10, #width
    file = dir,
    type = "jpeg", 
    bg = "transparent", 
    dpi = 300,
    units = "cm"
  )
   
  print(p)
  
  dev.off()
}






