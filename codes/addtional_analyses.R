
################################################################################
# additional analyses
################################################################################

# Result based on V16 and V23(stk.subset)
# V23, 10 subset splits

setwd("../../")
source("./developped cox/developped version 16, 5.08/cox.fun.R")
source("./functions/drfs.R")
source("./functions/drpca.R")
source("./functions/stk.cox_v4.R")
source("./functions/stk.glm_v3.R")
source("./functions/ipflasso.R")

library(glmnet)
library(survival)
library(dplyr)
library(tidyverse)

library(doParallel) 
library(doRNG) 
library(doSNOW) 

# all data
# datNames <- gsub(".RData", "", dir("../data/"))
# datNames <- datNames[-3]
# Parameters

paraList <- list(nsim   = 5,
                 nfolds = 10,
                 nCores = 25,
                 datNames = c("BLCA", "BRCA", "HNSC", "KIRC", "LAML", 
                              "LGG",  "LIHC", "LUAD", "PAAD", "SKCM"))
nsim   <- paraList$nsim
nfolds <- paraList$nfolds
nCores <- paraList$nCores
datNames <- paraList$datNames

#-------------------------------------------------------------------------------
# Explore the PCA probabilities
# Data, 260624
# update Lp of PCR method
#-------------------------------------------------------------------------------

#- predict.glmnet cannot handle the case: all coefficients are zero
getLp <- function(X, beta){
    
    beta <- drop(beta)
    
    if (length(beta) == ncol(X)) {
        # e.g., pcr
        lp <- X %*% beta
    } else {
        s = names(beta)
        if(length(s) > 0){
            lp <- X[,s, drop = F] %*% beta
        } else {
            lp <- matrix(0, nrow = nrow(X), ncol = 1)
        }
    }
    return(lp)
}


# set.seed(321)
for(nam in datNames) {
    
    print(nam)
    getData(nam, log2 = T, toBinary = F)
    
    n = nrow(X)
    p = ncol(X)
    
    #- get random splits, foldid, times
    getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
    time   <- y[ ,1]
    status <- y[ ,2]
    times <- quantile(time[status == 1], probs = seq(0.1, 0.9, length = 10))
    # timesRange <- quantile(time[status == 1], probs = c(0.1, 0.9))
    # times <- seq(timesRange[1], timesRange[2], length = 20)
    
    for(i in 1 : (nfolds * nsim)){
        
        samID <- ceiling(i / nfolds)
        # j, which fold of the i.nsim split being the test
        j <- i - (samID - 1) * nfolds
        # foldid for i.nsim splits and j fold
        sam <- sampleSplits[[samID]]
        load(paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV16.RData"))
        
        X.v <- X[sam == j,]
        y.v <- y[sam == j]
        
        X.t <- X[sam != j,]
        y.t <- y[sam != j]
        
        # mod PCR
        mod$pcr.las$lp <- getLp(X.t, mod$pcr.las$beta)
        save(list = "mod", 
             file = paste0("../output/", nam, "-", samID, "-", j, "-coxOptimV16.RData"))
        
    }
}


#-------------------------------------------------------------------------------
#  compare equal interval and quantile intervals
#-------------------------------------------------------------------------------



# set.seed(321)
for(nam in datNames) {
    
    print(nam)
    getData(nam, log2 = T, toBinary = F)
    
    n = nrow(X)
    p = ncol(X)
    
    #- get random splits, foldid, times
    getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
    time   <- y[ ,1]
    status <- y[ ,2]
    
    times <- quantile(time[status == 1], probs = seq(0.1, 0.9, length = 10))
    
    jpeg(paste("./", nam, ".jpeg"), width = 7, height = 3, units = "in", res = 300)
    
    par(mfrow = c(1, 2), mar = c(3.5, 3, 1.5, 1))
    plot(y, col = "black")
    abline(v = times2,col = "black", lty = "dashed")
    mtext("A, ten evenly distributed times", side = 3, line = 0.2, adj = 0)
    
    plot(y, col = "black")
    times2 <- seq(0, max(time), length = 10)
    abline(v = times, col = "black", lty = "dashed")
    mtext("B, ten quantiles of event times" , side = 3, line = 0.2, adj = 0)
    dev.off()

}




























