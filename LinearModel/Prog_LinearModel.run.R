## Packages

source("/code/LinearModel/Prog_LinearModel.R")

dir.create("/results/LM")
dir.create("/results/LM/Breast")
dir.create("/results/LM/Pan")
#############################################################
## linear model: breast cancer cell line
############################################################
load("/data/breast/DrugMeta.breast.zfixsample.RData")
load("/data/breast/DrugMeta.breast.fixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

z.aac <- list(mod.aac.ccle, mod.aac.ctrpv, mod.aac.gcsi,  mod.aac.gray, mod.aac.gdsc.v1, mod.aac.gdsc.v2, mod.aac.uhn)

z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data, z.exprs.gray.data, z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data, z.exprs.uhn.data)

study <- c("ccle", "ctrpv", "gcsi", "gray", "gdsc1", "gdsc2", "uhn")

for(i in 1:length(study)){
  
  print(i)
  lm.res <- lm.fun(exprs = z.exprs[[i]], aac= z.aac[[i]], drug = drugs)   
  save(lm.res, file = paste("/results/LM/Breast", 
                                     paste("lm.breast", study[i], "RData", sep="."), sep="/") ) 
  
}

#############################################################
## linear model: pan cancer cell line
############################################################
load("/data/pan/DrugMeta.pan.fixsample.mod.RData")
load("/data/pan/DrugMeta.pan.zfixsample.RData")

drugs <-c("Erlotinib","Lapatinib","Paclitaxel")

z.aac <- list(mod.aac.ccle, mod.aac.ctrpv, mod.aac.gcsi,
              mod.aac.gdsc.v1, mod.aac.gdsc.v2)

z.exprs <- list(z.exprs.ccle.data, z.exprs.ctrpv.data, z.exprs.gcsi.data, 
                z.exprs.gdsc.v1.data, z.exprs.gdsc.v2.data)

ccl.dat <- list(pan.ccl.ccle, pan.ccl.ctrpv, pan.ccl.gcsi, pan.ccl.gdsc.v1, pan.ccl.gdsc.v2)

study <- c("ccle", "ctrpv", "gcsi", "gdsc1", "gdsc2")

for(i in 1:length(study)){
  
  print(i)
  lm.res <- lm.adjusted.fun(exprs = z.exprs[[i]], aac= z.aac[[i]], drug = drugs, ccl = ccl.dat[[i]])   
  save(lm.res, file = paste("/results/LM/Pan", 
                            paste("lm.pan", study[i], "RData", sep="."), sep="/") ) 
  
}
