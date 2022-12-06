###############################################
## Packages
###############################################

library(meta)
library(GGally)
library(ggplot2)
library(metafor)
library(bayesmeta)
library(Cairo)
library(MASS)

dir.create("/results/simulation")
###############################################
## parameters or factors
###############################################

method <- c("DL")
n.study <- c(3, 7, 10)
beta0 <-  c( 0.2, 0.5, 0.8)
tau.val <- round(sqrt(c(0.001,0.20,0.5, 0.80)),2) 
sigma.val <- sqrt(c( 0.1, 0.45, 1)) 
# Set the number of simulation to 1000
n.sim <- 50 
rao.val <- c(0.1, 0.7) ## for one duplication

#############################################################################
## simulate data: one duplication
#############################################################################
# repeat the following simulation for one, two and three duplications
set.seed(135)

meta.sim.data <- lapply(1:n.sim, function(l){
  
  
  sim.data <- lapply(1:length(n.study), function(i){
    
    
    sim.sigma.data <- lapply(1:length(sigma.val), function(j){
      
      
      sim.tau.data <- lapply(1:length(tau.val), function(k){
      
        print(k)
        
        res.rao <- lapply(1:length(rao.val), function(s){
        
          print(s)
          
        tau.matrix <- diag(rep(tau.val[k]^2, n.study[i]))
        
        ## one duplication
        tau.matrix[1,2] <- rao.val[s] * sigma.val[j] * sigma.val[j]
        tau.matrix[2,1] <- rao.val[s] * sigma.val[j] * sigma.val[j]
        
        
        sigma.matrix <- diag(rep(sigma.val[j]^2, n.study[i]))
      
        res.method <- lapply(1:length(beta0), function(m){
          
        y.obs <- beta0[m] + mvrnorm(n=1, mu = rep(0, n.study[i]), Sigma = tau.matrix + sigma.matrix)
        
        se.obs <- rep(sigma.val[j], n.study[i])
          
          print(m)
          res <- metagen(y.obs, se.obs, method.tau = "DL")
          
          cov.ci <- ifelse(res$lower.random < beta0[m] & res$upper.random > beta0[m], 1, 0)
          res.val <- data.frame(sim.id = l,
                                method = "DL", 
                                n.study = n.study[i],
                                tau0 = tau.val[k],
                                rao = rao.val[s],
                                beta0 = beta0[m],
                                sigma.val = sigma.val[j],
                                TE.random = res$TE.random, 
                                lower.random = res$lower.random, upper.random = res$upper.random,
                                length.ci = res$upper.random - res$lower.random,
                                tau = res$tau,
                                cov.ci = cov.ci)
          
          res.val
          
        })
          
        do.call(rbind, res.method)
        
      
        })
        
        do.call(rbind, res.rao)
        
      })
      
      do.call(rbind, sim.tau.data)
      
    })
    
    do.call(rbind, sim.sigma.data)
    
  })
  
  do.call(rbind, sim.data)
  
})

meta.sim.data <- do.call(rbind, meta.sim.data)

save(meta.sim.data, 
     file = "/results/simulation/res.dept.one.duplication.RData")

#############################################################################
## simulate data: two duplication
#############################################################################
set.seed(135)

meta.sim.data <- lapply(1:n.sim, function(l){
  
  
  sim.data <- lapply(1:length(n.study), function(i){
    
    
    sim.sigma.data <- lapply(1:length(sigma.val), function(j){
      
      
      sim.tau.data <- lapply(1:length(tau.val), function(k){
        
        
        res.rao <- lapply(1:length(rao.val), function(s){
          
          print(s)
          
          tau.matrix <- diag(rep(tau.val[k]^2, n.study[i]))
          
          ## one duplication
          tau.matrix[1,2] <- rao.val[s] * sigma.val[j] * sigma.val[j]
          tau.matrix[2,1] <- rao.val[s] * sigma.val[j] * sigma.val[j]
          
          ## two duplications using rao = 0.2
          
          tau.matrix[2,3] <- 0.2 * sigma.val[j] * sigma.val[j]
          tau.matrix[3,2] <- 0.2 * sigma.val[j] * sigma.val[j]

          
          sigma.matrix <- diag(rep(sigma.val[j]^2, n.study[i]))
          
          res.method <- lapply(1:length(beta0), function(m){
            
            y.obs <- beta0[m] + mvrnorm(n=1, mu = rep(0, n.study[i]), Sigma = tau.matrix + sigma.matrix)
            
            se.obs <- rep(sigma.val[j], n.study[i])
            
            print(m)
            res <- metagen(y.obs, se.obs, method.tau = "DL")
            
            cov.ci <- ifelse(res$lower.random < beta0[m] & res$upper.random > beta0[m], 1, 0)
            res.val <- data.frame(sim.id = l,
                                  method = "DL", 
                                  n.study = n.study[i],
                                  tau0 = tau.val[k],
                                  rao = rao.val[s],
                                  beta0 = beta0[m],
                                  sigma.val = sigma.val[j],
                                  TE.random = res$TE.random, 
                                  lower.random = res$lower.random, upper.random = res$upper.random,
                                  length.ci = res$upper.random - res$lower.random,
                                  tau = res$tau,
                                  cov.ci = cov.ci)
            
            res.val
            
          })
          
          do.call(rbind, res.method)
          
          
        })
        
        do.call(rbind, res.rao)
        
      })
      
      do.call(rbind, sim.tau.data)
      
    })
    
    do.call(rbind, sim.sigma.data)
    
  })
  
  do.call(rbind, sim.data)
  
})

meta.sim.data <- do.call(rbind, meta.sim.data)

save(meta.sim.data, 
     file = "/results/simulation/res.dept.two.duplication.RData")


#############################################################################
## simulate data: three duplication
#############################################################################
set.seed(135)

meta.sim.data <- lapply(1:n.sim, function(l){
  
  sim.data <- lapply(1:length(n.study), function(i){
    
    sim.sigma.data <- lapply(1:length(sigma.val), function(j){
      
      sim.tau.data <- lapply(1:length(tau.val), function(k){
        
        res.rao <- lapply(1:length(rao.val), function(s){
          
          tau.matrix <- diag(rep(tau.val[k]^2, n.study[i]))
          
          ## one duplication
          tau.matrix[1,2] <- rao.val[s] * sigma.val[j] * sigma.val[j]
          tau.matrix[2,1] <- rao.val[s] * sigma.val[j] * sigma.val[j]
          
          ## two duplications using rao = 0.2
          
          tau.matrix[2,3] <- 0.2 * sigma.val[j] * sigma.val[j]
          tau.matrix[3,2] <- 0.2 * sigma.val[j] * sigma.val[j]
          
          ## three duplications using rao = 0.6
          
          tau.matrix[1,3] <- 0.6 * sigma.val[j] * sigma.val[j]
          tau.matrix[3,1] <- 0.6 * sigma.val[j] * sigma.val[j]
          
          sigma.matrix <- diag(rep(sigma.val[j]^2, n.study[i]))
          
          res.method <- lapply(1:length(beta0), function(m){
            
            y.obs <- beta0[m] + mvrnorm(n=1, mu = rep(0, n.study[i]), Sigma = tau.matrix + sigma.matrix)
            
            se.obs <- rep(sigma.val[j], n.study[i])
            
            print(m)
            res <- metagen(y.obs, se.obs, method.tau = "DL")
            
            cov.ci <- ifelse(res$lower.random < beta0[m] & res$upper.random > beta0[m], 1, 0)
            res.val <- data.frame(sim.id = l,
                                  method = "DL", 
                                  n.study = n.study[i],
                                  tau0 = tau.val[k],
                                  rao = rao.val[s],
                                  beta0 = beta0[m],
                                  sigma.val = sigma.val[j],
                                  TE.random = res$TE.random, 
                                  lower.random = res$lower.random, upper.random = res$upper.random,
                                  length.ci = res$upper.random - res$lower.random,
                                  tau = res$tau,
                                  cov.ci = cov.ci)
            
            res.val
            
          })
          
          do.call(rbind, res.method)
          
          
        })
        
        do.call(rbind, res.rao)
        
      })
      
      do.call(rbind, sim.tau.data)
      
    })
    
    do.call(rbind, sim.sigma.data)
    
  })
  
  do.call(rbind, sim.data)
  
})

meta.sim.data <- do.call(rbind, meta.sim.data)

save(meta.sim.data, 
     file = "/results/simulation/res.dept.three.duplication.RData")

##################################################################
##################################################################
## Combine all three duplications results
##################################################################
##################################################################

load("/results/simulation/res.dept.one.duplication.RData") 

method <- unique(meta.sim.data $method)
n.study <- unique(meta.sim.data$n.study)
tau.val <- unique(meta.sim.data$tau0)
sigma.val <- unique(meta.sim.data$sigma.val)
beta0 <- unique(meta.sim.data$beta0)
rao.val <- unique(meta.sim.data$rao)

bias.mse.res <- lapply(1:length(beta0), function(m){
  
  rao.res <- lapply(1:length(rao.val), function(s){
    
    method.res <-  lapply(1:length(n.study), function(k){
      
      tau.res <- lapply(1:length(tau.val), function(i){
        
        sigma.res <- lapply(1:length(sigma.val), function(j){
          
          res <- meta.sim.data[meta.sim.data$n.study == n.study[k] &
                                 meta.sim.data$tau0 == tau.val[i] &
                                 meta.sim.data$sigma.val == sigma.val[j] &
                                 meta.sim.data$rao == rao.val[s] &
                                 meta.sim.data$beta0 == beta0[m],  ]
          
          sd.random <- sd(res$TE.random)
          sd.tau <- sd(res$tau)
          mean.random <- mean(res$TE.random)
          mean.tau <- mean(res$tau)
          bias.random <- (mean(res$TE.random) - res$beta0)
          mse.random <- bias.random^2 + (sd.random)^2
          bias.tau <- (mean(res$tau) - res$tau0)^2
          mse.tau <- bias.tau^2 + (sd.tau)^2
          sd.ci <- sd(res$length.ci)
          mean.ci <- mean(res$length.ci)
          r.bias.random <- mean(res$TE.random)/res$beta0
          r.mse.random <- mse.random/(res$beta0)^2
          r.bias.tau <- mean(res$tau)/res$tau0
          r.mse.tau <- mse.tau/(res$tau0)^2
          cov.rate <- nrow(res[res$cov.ci == 1, ])/nrow(res)
          
          data.frame(n.study = n.study[k], 
                     beta0 = beta0[m],
                     rao = rao.val[s],
                     tau0 = tau.val[i],
                     sigma.val = sigma.val[j],
                     bias.random,
                     mse.random,
                     bias.tau,
                     mse.tau,
                     sd.random,
                     sd.tau,
                     mean.random,
                     mean.tau,
                     sd.ci ,
                     mean.ci,
                     r.bias.random = r.bias.random,
                     r.mse.random = r.mse.random,
                     r.bias.tau = r.bias.tau,
                     r.mse.tau = r.mse.tau,
                     cov.rate)
          
        })
        
        do.call(rbind, sigma.res)
        
      })
      
      do.call(rbind, tau.res)
      
    })
    
    do.call(rbind, method.res)
    
  })
  
  do.call(rbind, rao.res)
  
})
  
bias.mse.res <- do.call(rbind, bias.mse.res)
bias.mse.res$tau0 <- round(bias.mse.res$tau0,2)
bias.mse.res$tau0  <- factor(bias.mse.res$tau0)
bias.mse.res$sigma.val <- round(bias.mse.res$sigma.val, 2)
bias.mse.res$sigma.val <- factor(bias.mse.res$sigma.val)
bias.mse.res$n.study <- factor(bias.mse.res$n.study)

bias.mse.res$rao <- factor(bias.mse.res$rao)
bias.mse.res.one <- data.frame(duplication = "one", bias.mse.res)


load("/results/simulation/res.dept.two.duplication.RData") 

method <- unique(meta.sim.data $method)
n.study <- unique(meta.sim.data$n.study)
tau.val <- unique(meta.sim.data$tau0)
sigma.val <- unique(meta.sim.data$sigma.val)
beta0 <- unique(meta.sim.data$beta0)
rao.val <- unique(meta.sim.data$rao)

bias.mse.res <- lapply(1:length(beta0), function(m){
  
  rao.res <- lapply(1:length(rao.val), function(s){
    
    method.res <-  lapply(1:length(n.study), function(k){
      
      tau.res <- lapply(1:length(tau.val), function(i){
        
        sigma.res <- lapply(1:length(sigma.val), function(j){
          
          res <- meta.sim.data[meta.sim.data$n.study == n.study[k] &
                                 meta.sim.data$tau0 == tau.val[i] &
                                 meta.sim.data$sigma.val == sigma.val[j] &
                                 meta.sim.data$rao == rao.val[s] &
                                 meta.sim.data$beta0 == beta0[m],  ]
          
          sd.random <- sd(res$TE.random)
          sd.tau <- sd(res$tau)
          mean.random <- mean(res$TE.random)
          mean.tau <- mean(res$tau)
          bias.random <- (mean(res$TE.random) - res$beta0)
          mse.random <- bias.random^2 + (sd.random)^2
          bias.tau <- (mean(res$tau) - res$tau0)^2
          mse.tau <- bias.tau^2 + (sd.tau)^2
          sd.ci <- sd(res$length.ci)
          mean.ci <- mean(res$length.ci)
          r.bias.random <- mean(res$TE.random)/res$beta0
          r.mse.random <- mse.random/(res$beta0)^2
          r.bias.tau <- mean(res$tau)/res$tau0
          r.mse.tau <- mse.tau/(res$tau0)^2
          cov.rate <- nrow(res[res$cov.ci == 1, ])/nrow(res)
          
          data.frame(n.study = n.study[k], 
                     beta0 = beta0[m],
                     rao = rao.val[s],
                     tau0 = tau.val[i],
                     sigma.val = sigma.val[j],
                     bias.random,
                     mse.random,
                     bias.tau,
                     mse.tau,
                     sd.random,
                     sd.tau,
                     mean.random,
                     mean.tau,
                     sd.ci ,
                     mean.ci,
                     r.bias.random = r.bias.random,
                     r.mse.random = r.mse.random,
                     r.bias.tau = r.bias.tau,
                     r.mse.tau = r.mse.tau,
                     cov.rate)
          
        })
        
        do.call(rbind, sigma.res)
        
      })
      
      do.call(rbind, tau.res)
      
    })
    
    do.call(rbind, method.res)
    
  })
  
  do.call(rbind, rao.res)
  
})

bias.mse.res <- do.call(rbind, bias.mse.res)
bias.mse.res$tau0 <- round(bias.mse.res$tau0,2)
bias.mse.res$tau0  <- factor(bias.mse.res$tau0)
bias.mse.res$sigma.val <- round(bias.mse.res$sigma.val, 2)
bias.mse.res$sigma.val <- factor(bias.mse.res$sigma.val)
bias.mse.res$n.study <- factor(bias.mse.res$n.study)

bias.mse.res$rao <- factor(bias.mse.res$rao)
bias.mse.res.two <- data.frame(duplication = "two", bias.mse.res)


load("/results/simulation/res.dept.three.duplication.RData") 

method <- unique(meta.sim.data $method)
n.study <- unique(meta.sim.data$n.study)
tau.val <- unique(meta.sim.data$tau0)
sigma.val <- unique(meta.sim.data$sigma.val)
beta0 <- unique(meta.sim.data$beta0)
rao.val <- unique(meta.sim.data$rao)

bias.mse.res <- lapply(1:length(beta0), function(m){
  
  rao.res <- lapply(1:length(rao.val), function(s){
    
    method.res <-  lapply(1:length(n.study), function(k){
      
      tau.res <- lapply(1:length(tau.val), function(i){
        
        sigma.res <- lapply(1:length(sigma.val), function(j){
          
          res <- meta.sim.data[meta.sim.data$n.study == n.study[k] &
                                 meta.sim.data$tau0 == tau.val[i] &
                                 meta.sim.data$sigma.val == sigma.val[j] &
                                 meta.sim.data$rao == rao.val[s] &
                                 meta.sim.data$beta0 == beta0[m],  ]
          
          sd.random <- sd(res$TE.random)
          sd.tau <- sd(res$tau)
          mean.random <- mean(res$TE.random)
          mean.tau <- mean(res$tau)
          bias.random <- (mean(res$TE.random) - res$beta0)
          mse.random <- bias.random^2 + (sd.random)^2
          bias.tau <- (mean(res$tau) - res$tau0)^2
          mse.tau <- bias.tau^2 + (sd.tau)^2
          sd.ci <- sd(res$length.ci)
          mean.ci <- mean(res$length.ci)
          r.bias.random <- mean(res$TE.random)/res$beta0
          r.mse.random <- mse.random/(res$beta0)^2
          r.bias.tau <- mean(res$tau)/res$tau0
          r.mse.tau <- mse.tau/(res$tau0)^2
          cov.rate <- nrow(res[res$cov.ci == 1, ])/nrow(res)
          
          data.frame(n.study = n.study[k], 
                     beta0 = beta0[m],
                     rao = rao.val[s],
                     tau0 = tau.val[i],
                     sigma.val = sigma.val[j],
                     bias.random,
                     mse.random,
                     bias.tau,
                     mse.tau,
                     sd.random,
                     sd.tau,
                     mean.random,
                     mean.tau,
                     sd.ci ,
                     mean.ci,
                     r.bias.random = r.bias.random,
                     r.mse.random = r.mse.random,
                     r.bias.tau = r.bias.tau,
                     r.mse.tau = r.mse.tau,
                     cov.rate)
          
        })
        
        do.call(rbind, sigma.res)
        
      })
      
      do.call(rbind, tau.res)
      
    })
    
    do.call(rbind, method.res)
    
  })
  
  do.call(rbind, rao.res)
  
})

bias.mse.res <- do.call(rbind, bias.mse.res)
bias.mse.res$tau0 <- round(bias.mse.res$tau0,2)
bias.mse.res$tau0  <- factor(bias.mse.res$tau0)
bias.mse.res$sigma.val <- round(bias.mse.res$sigma.val, 2)
bias.mse.res$sigma.val <- factor(bias.mse.res$sigma.val)
bias.mse.res$n.study <- factor(bias.mse.res$n.study)

bias.mse.res$rao <- factor(bias.mse.res$rao)
bias.mse.res.three <- data.frame(duplication = "three", bias.mse.res)


bias.mse.res <- rbind(bias.mse.res.one, bias.mse.res.two, bias.mse.res.three)

save(bias.mse.res, 
     file = "/results/simulation/res.dept.duplication.all.RData")


##############################################################################
## Visualization
##############################################################################
load("/results/simulation/res.dept.duplication.all.RData")

beta0 <- c(0.2, 0.5, 0.8)
rao <- c(0.1, 0.7)   

bias.mse.res <- bias.mse.res[bias.mse.res$sigma.val != 0.1, ]

for(i in 1:length(beta0)){
  
  for(j in 1:length(rao)){
    
    dat <- bias.mse.res[bias.mse.res$beta0 == beta0[i] & # 0.2, 0.5, 0.8
                          bias.mse.res$rao == rao[j], # 0.1, 0.7 
                        ]
    
    dir <- paste("/results/simulation/", paste("mse", paste("beta", beta0[i], sep=""), 
                                                      paste("rao", rao[j], sep=""), sep = "_"), ".jpeg", 
                 sep="")
    
    Cairo::Cairo(
      20, #length
      15, #width
      file = dir,
      type = "jpeg", 
      bg = "transparent", 
      dpi = 300,
      units = "cm"
    )
    
    
    p <- ggplot(dat, 
           aes(x= n.study , y = r.mse.random,   # mse.random
               col = duplication, group=duplication)) + 
      geom_point() +
      geom_line() +
      facet_grid(sigma.val ~ tau0 ) +
      scale_colour_manual(values=c("#542788", "#9e0142", "#66c2a5")) + 
      ylab(expression(paste("relative MSE of true effect ", "(",beta, ")"))) +
      ggtitle(expression(paste(beta, "=", 0.8 ,", ", rho,'=', 0.7))) +
      xlab(" ") +
      theme(axis.text.x=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            axis.text.y=element_text(size=12),
            strip.text = element_text(size=14, face="bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            legend.position="bottom",
            legend.text = element_text(size = 12, face="bold"),
            legend.title = element_blank(),
            plot.title=element_text(hjust=0.5))
    
    print(p)
    dev.off()
    
    dir <- paste("/results/simulation/", paste("cp", paste("beta", beta0[i], sep=""), 
                                                      paste("rao", rao[j], sep=""), sep = "_"), ".jpeg", 
                 sep="")
    
    Cairo::Cairo(
      20, #length
      15, #width
      file = dir,
      type = "jpeg", 
      bg = "transparent", 
      dpi = 300,
      units = "cm"
    )
    
    p <- ggplot(dat, 
           aes(x= n.study , y = cov.rate,   # mse.random
               col = duplication, group=duplication)) + 
      geom_point() +
      geom_line() +
      facet_grid(sigma.val ~ tau0 ) +
      scale_colour_manual(values=c("#542788", "#9e0142", "#66c2a5")) + 
      ylab("coverage probability") +
      ggtitle(expression(paste(beta, "=", 0.8 ,", ", rho,'=', 0.7))) +
      xlab(" ") +
      theme(axis.text.x=element_text(size=14),
            axis.title=element_text(size=14,face="bold"),
            axis.text.y=element_text(size=12),
            strip.text = element_text(size=14, face="bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            plot.background = element_blank(), 
            axis.line = element_line(colour = "black"),
            legend.position="bottom",
            legend.text = element_text(size = 10, face="bold"),
            legend.title = element_blank(),
            plot.title=element_text(hjust=0.5))
    
    print(p)
    
    dev.off()

  }
  
}
