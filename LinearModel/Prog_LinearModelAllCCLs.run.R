## Packages

source("/code/LinearModel/Prog_LinearModel.R")

dir.create("/results/LM/Breast/AllCCLs")
#############################################################
## linear model: breast cancer cell line
############################################################
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

