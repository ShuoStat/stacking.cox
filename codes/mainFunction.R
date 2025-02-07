
################################################################################
# based on version-29
################################################################################

# load packages and functions

source("./cox.fun.R")
source("./helper.R")
source("./drfs.R")
source("./drpca.R")
source("./stk.cox_v7.R")
source("./stk.glm_v4.R")
source("./ipflasso.R")

# load package

packages <- c("glmnet", "survival", "dplyr", "tidyverse", "doParallel", "doRNG", "doSNOW")
loadPackages(packages)

# parameters

paraList <- list(nsim   = 10,
                 nfolds = 10,
                 nCores = 25, 
                 datNames = c("BLCA", "BRCA", "HNSC", "KIRC", "LAML", 
                              "LGG",  "LIHC", "LUAD", "PAAD", "SKCM"))
nsim   <- paraList$nsim # Number of repetitions for k-fold cross-validation
nfolds <- paraList$nfolds # Number of folds in k-fold cross-validation
nCores <- paraList$nCores # Number of CPU cores allocated for parallel processing (e.g., 25 cores)
datNames <- paraList$datNames # List of dataset names to be used

#-------------------------------------------------------------------------------
#- Get clinical information
#-------------------------------------------------------------------------------

sumDat <- c()

for(nam in datNames){

  cat(nam, ":", as.character(Sys.time()), "\n\n")
  
  # get data, output X, y, and blocks to global environment
  getData(nam, log2 = T, toBinary = F)
  print(sum(y[,2]))
  
  # get N clin
  N_clin <- sum(grepl("^clin_", colnames(X)))
  N_mol  <- sum(grepl("^m", colnames(X)))
  
  # get sampel size
  N <- nrow(X)
  
  # get N events
  N_event <- sum(y[,2])
  
  # get median survival time
  surv_fit <- survfit(y ~ 1)
  Meidan_S <- round(quantile(surv_fit$time, 0.5))
  
  inF <- setNames(c(nam, N, N_event, N_clin, N_mol, Meidan_S),
                  c("Data", "Sample Size", "Event", "N clinical variables", "N molecular variables", "Median Survival"))

  sumDat <- rbind(sumDat, inF)
}

print(sumDat)
write.csv(sumDat, file ="../results/basicInfo.csv")

#-------------------------------------------------------------------------------
# run it: training models based on 10 * 10-fold-cross-validation
#-------------------------------------------------------------------------------

# run the selected models
# basic models: clinical, molecular, naive lasso, dimension reduction based on feature screening and pca, ipflasso
# stk, stk.pf and stk.cp
# stkPF.pro, subset stacking
whichRun <- c("basic", "stk", "stkPF.pro")

for(nam in datNames){
  
  # obtain data
  cat(nam, ":", as.character(Sys.time()), "\n\n")
  getData(nam, log2 = T, toBinary = F)
  
  # sample size and feature size
  n <- nrow(X)
  p <- ncol(X)
  
  # get split for the data
  getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
  time   <- y[,1]
  status <- y[,2]
  
  times <- quantile(time[status == 1], probs = seq(0.1, 0.9, length = 10))
  cox.time <- quantile(time[status == 1], probs = 0.5)
  
  # using 500 features to examine the codes
  # X <- X[,1:500]
  # blocks <- blocks[1:500]
  
  # parallel runing setup
  Cores <- pmin(nCores, detectCores() - 1)
  cl <- makeCluster(Cores)
  registerDoSNOW(cl)
  
  # initialize progress bar
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  # perform parallel cross-validation using foreach
  foreach(i = 1 : (nfolds * nsim),
          .packages = c("glmnet", "survival"),
          .options.snow = opts,
          .verbose = TRUE) %dopar% {
            
            # determine which simulation and fold is currently running
            samID <- ceiling(i / nfolds) # current simulation ID
            j <- i - (samID - 1) * nfolds # current fold ID
            
            # retrieve train-test split information
            sam <- sampleSplits[[samID]]
            foldid <- foldids[[samID]][[j]]
            foldid.internal <- foldid.internals[[samID]][[j]]
            
            # define training and validation sets
            X.t <- X[sam != j,]  # training features
            y.t <- y[sam != j]    # training outcome
            X.v <- X[sam == j,]  # validation features
            y.v <- y[sam == j]    # validation outcome

            # argument list for model training
            argList <- list(X = X.t, y = y.t, 
                            times = times,
                            cox.time = cox.time, 
                            blocks = blocks,
                            foldid = foldid,
                            intercept = TRUE,
                            foldid.internal = foldid.internal)
            
            
            # train and save the basic model if selected
            if ("basic" %in% whichRun) {
              ncomp <- nrow(X.t) - 1
              argList$ncomp <- ncomp
              mod <- do.call(fit.basic, argList)
              save(list = "mod", file = paste0("../output/basic-", nam, "-", samID, "-", j, ".RData"))
            }
            
            # train and save the stacking model if selected
            if("stk" %in% whichRun){
              mod <- do.call(fit.stk, argList)
              save(list = "mod", file = paste0("../output/stk-", nam, "-", samID, "-", j, ".RData"))
            }
            
            # train and save the subset stacking if selected
            if ("stkPF.pro" %in% whichRun) {
              mod <- do.call(fit.stkPF.pro, argList)
              save(list = "mod", file = paste0("../output/stkPF.pro-", nam, "-", samID, "-", j, ".RData"))
            }
          }
  
  stopCluster(cl)
}


#-------------------------------------------------------------------------------
#- prediction 
#-------------------------------------------------------------------------------

# indicate whether sub-model predictions use averaged cross-validation sub-models. Here is TRUE
ifCVMods = TRUE

ibss <- list() # list for ibss results 
aucs <- list() # list for aucs results 
logLiks <- list() # list for logLik 

# set.seed(321) 
for(nam in datNames) { 
  
  print(nam) 
  getData(nam, log2 = T, toBinary = F) 
  
  n = nrow(X) 
  p = ncol(X) 
  
  # # try
  # X <- X[,1:500]
  # blocks <- blocks[1:500]

  #- get random splits, foldid, times
  getSplits(n, nsim, nfolds, seed = 1248372, foldid.internal = T)
  time   <- y[ ,1]
  status <- y[ ,2]
  
  # select 20 time points evenly spaced between the 0.1 and 0.9 quantiles of observed event times
  timesRange <- quantile(time[status == 1], probs = c(0.1, 0.9))
  times <- seq(timesRange[1], timesRange[2], length = 20)
  
  # parallel runing
  Cores <- pmin(nCores, detectCores() - 1)
  cl <- makeCluster(Cores)
  registerDoSNOW(cl)
  
  # initialize progress bar
  pb <- txtProgressBar(min = 1, max = (nfolds * nsim), style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)
  
  res <- foreach(i = 1 : (nfolds * nsim), 
                 .packages = c("dplyr", "glmnet", "survival"),
                 .options.snow = opts,
                 .verbose = TRUE) %dopar% {
                   
                   samID <- ceiling(i / nfolds)
                   # j, which fold of the i.nsim split being the test
                   j <- i - (samID - 1) * nfolds
                   # foldid for i.nsim splits and j fold
                   sam <- sampleSplits[[samID]]
                   
                   load(paste0("../output/basic_", nam, "_", samID, "_", j, ".RData")) 
                   mods <- mod  
                   load(paste0("../output/stk_", nam, "_", samID, "_", j, ".RData"))
                   mods <- c(mods, mod[-1]) # mod[1] is clinical model for checking the identity between models
                   load(paste0("../output/stkPF.pro_", nam, "_", samID, "_", j, ".RData"))
                   mod <- c(mods, mod[-1]) 

                   # if a single model return (e.g., stkPF.pro), transform it into a list. 
                   # the format used for the following analysis is a list containing models
                   if (length(mod) == 1) 
                     mod = list(mod) 
                   
                   # training and validation data
                   X.v <- X[sam == j,] 
                   y.v <- y[sam == j]
                   
                   X.t <- X[sam != j,]
                   y.t <- y[sam != j]
                   # derive the prediction probabilities
                   probs <- predProbs(X.v, mod, times, ifCVMods = ifCVMods)

                   # flatten list
                   flattenList <- function(obj, target) {
                     new <- obj[[target]]
                     names(new) <- paste0(target, ":", names(new))
                     obj <- append(obj, new)
                     obj[[target]] <- NULL
                     obj
                   }
                   
                   for(k in names(probs)[grep("stk", names(probs))]) {
                     probs <- flattenList(probs, target = k)
                   }
                   
                   # get Time 
                   ibsTime <- lapply(probs, getTimeIBS,
                                     Xclin.t = X.t[ ,blocks == 1],
                                     Xclin.v = X.v[ ,blocks == 1],
                                     y.t = y.t,
                                     y.v = y.v,
                                     times = times,
                                     cens.method = "clinicalAdjusted",
                                     truncICPW = 0.95,
                                     sumRes = FALSE)
                   
                   # get AUC
                   aucTime <- lapply(probs, getTimeAUC,
                                     Xclin.t = X.t[ ,blocks == 1],
                                     Xclin.v = X.v[ ,blocks == 1],
                                     y.t = y.t,
                                     y.v = y.v,
                                     times = times,
                                     cens.method = "clinicalAdjusted",
                                     truncICPW = 0.95,
                                     sumRes = FALSE)
                   
                   # get LogLik
                   logLikTime <- lapply(probs, getTimelogLik,
                                        y.v = y.v,
                                        times = times,
                                        sumRes = FALSE)
                   
                   return(list(ibsTime = ibsTime,
                               aucTime = aucTime,
                               logLikTime = logLikTime))
                 }
  # from res to ibss, aucs, logLiks
  ibss[[nam]] <- lapply(res, `[[`, "ibsTime")
  aucs[[nam]] <- lapply(res, `[[`, "aucTime")
  logLiks[[nam]] <- lapply(res, `[[`, "logLikTime")
  
  stopCluster(cl)
}

save(list = c("aucs", "ibss", "logLiks"), file = paste0("../output/pred.RData"))

#-------------------------------------------------------------------------------
# Summarize the prediction 
#-------------------------------------------------------------------------------

# function to summarize results
sumFun <- function(objs, fun = colMeans) {
  
  # objs, aucs, ibss, or logLiks
  areaSum <- function(obj) {
    obj <- na.omit(obj)
    tRange <- abs(range(obj[ ,"time"])[1] - range(obj[ ,"time"])[2])
    trapezoidal(obj[,"time"], obj[,2]) / tRange
    # sum(obj[,2])
  }
  
  # function to summarize results for each object
  sumRes <- function(obj, fun){
    res <- lapply(obj, function(x) unlist(lapply(x, areaSum)))
    res <- Reduce(rbind, res)
    do.call(fun, list(res))
  }
  
  tab <- lapply(objs, sumRes, fun = fun)
  tab <- Reduce(rbind, tab)
  rownames(tab) <- names(objs)
  return(tab)
}

# summarize results
sumResults <- function(target, ifSE = FALSE){
  
  # target, which to run, between aucs, ibss, and logLiks
  load(paste0("../output/pred.RData"))
  
  # summarize the data
  reMat <- sumFun(get(target), fun = colMeans)
  
  # how many replicates
  n <- length(get(target)$BLCA)
  
  if (ifSE) {
    reMatSe <- sumFun(get(target), fun = function(x) apply(x, 2, sd) / sqrt(n))
  }
  
  # calculate mean
  re <- as.data.frame(reMat) %>% 
    select(clin, mol, naive, fs.las, pcr.las, ipf.las, `stk.pf:ibsLoss`, `stk.pf:logLik:ridge`, `stk.pf:PaLogLik:ridge`, `stk.comp:ibsLoss`,
           `stk.comp:logLik:ridge`, `stk.comp:PaLogLik:ridge`) %>%
    rename_all( ~ c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPFLasso", "Stk(PF,IBS)", "Stk(PF,NLL)", "Stk(PF,PaNLL)",
                    "Stk(CP,IBS)", "Stk(CP,NLL)", "Stk(CP,PaNLL)")) %>%
    mutate(across(everything(), ~ formatC(.x, format = "f", digits = 3)))
  
  # calculate standard error
  if (ifSE) {
    re.se <- as.data.frame(reMatSe) %>% 
      select(clin, mol, naive, fs.las, pcr.las, ipf.las, `stk.pf:ibsLoss`, `stk.pf:logLik:ridge`, `stk.pf:PaLogLik:ridge`, `stk.comp:ibsLoss`,
             `stk.comp:logLik:ridge`, `stk.comp:PaLogLik:ridge`) %>%
      rename_all( ~ c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPFLasso", "Stk(PF,IBS)", "Stk(PF,NLL)", "Stk(PF,PaNLL)",
                      "Stk(CP,IBS)", "Stk(CP,NLL)", "Stk(CP,PaNLL)")) %>%
      mutate(across(everything(), ~ formatC(.x, format = "f", digits = 3)))
    
    re <- mapply(function(mean, se) {
      paste0(mean, "(", se, ")")
    }, mean = re, se = re.se)
    rownames(re) <- rownames(re.se)
  }
  
  # rowname data to first col
  re <- tibble::rownames_to_column(as.data.frame(re), var = "Data")
  openxlsx::write.xlsx(re, file = paste0("../results/", target, ".xlsx"), colNames = T)
  
  print(re)
  return(re)
}

# summarize the results
sumResults(target = "ibss", ifSE = TRUE)
sumResults(target = "aucs", ifSE = TRUE)
sumResults(target = "logLiks", ifSE = TRUE)

#-------------------------------------------------------------------------------
# Line Plots: AUC by Time
#-------------------------------------------------------------------------------

# function to generate line plots for AUC over time
Plot_Time <- function(obj, ylab, title) {
  
  # loop through each model's data and merge them into a single data.frame
  c_obj <- c()
  for(i in seq_along(obj)){
    tmp_obj <- obj[[i]]
    tmp_obj <- cbind(as.data.frame(Reduce(rbind, tmp_obj)), mod = rep(names(tmp_obj), each = nrow(tmp_obj[[1]])))
    c_obj   <- rbind(c_obj, tmp_obj)
  }
  
  colnames(c_obj) = c("time", "value", "mod")
  
  datPlot <- c_obj %>%
    filter(value != "NaN") %>%
    filter(mod %in% c("naive", "ipf.las", "stk.pf:logLik:ridge")) %>%
    mutate(mod = recode(mod, "naive" = "Lasso", "ipf.las" = "IPFLasso", "stk.pf:logLik:ridge" = "Stk(PF,NLL)"),
           mod = factor(mod, levels = c("Lasso", "IPFLasso", "Stk(PF,NLL)"))) %>%
    group_by(time, mod) %>%
    summarise(value = mean(value))
  
  ggplot2::ggplot(datPlot, aes(x = time, y = value, color  = mod)) +
    geom_point(size = 1) +
    geom_line() +
    theme_bw() +
    labs(title = title,
         x = "Time,Days",
         y = ylab,
         color = "Model") + 
    theme(legend.position = "bottom", 
          legend.title = element_blank())
}

# load prediction results
load(paste0("../output/pred.RData"))
# loop through dataset names and generate plots for each
plotList <- list()
for(i in datNames) {
  plotList[[i]] <- Plot_Time(aucs[[i]], ylab = "AUC", title = i)
}

# save plots
p <- ggpubr::ggarrange(plotlist = plotList, nrow = 2, ncol = 5, common.legend = TRUE, legend = "bottom")
ggsave(filename = "../results/timeAUC.pdf", width = 12, height = 5)

#-------------------------------------------------------------------------------
# Subset Stacking 
#-------------------------------------------------------------------------------

vsSS <- function(target, ifSE = FALSE, to = c("table", "figure")) {
  # to = c("table", "figure"), output tables and figures
  
  load(paste0("../output/pred.RData"))
  n <- length(get(target)[[1]])
  
  # remove duplicated clin
  reMat <- sumFun(get(target), fun = colMeans)
  reMatSe <- sumFun(get(target), fun = function(x) apply(x, 2, sd) / sqrt(n))
  
  # selected moehods for comparison
  # to table
  sel <- c("stk.pf:ibsLoss", "stk.subset:ibsLoss", "stk.pf:logLik:ridge", "stk.subset:logLik:ridge", "stk.pf:PaLogLik:ridge", "stk.subset:PaLogLik:ridge")
  reNames <- c("Stk(PF,IBS)", "Stk(SS,IBS)", "Stk(PF,NLL)", "Stk(SS,NLL)", "Stk(PF,PaNLL)", "Stk(SS,PaNLL)")
  # I cannot use select because of `:`
  re <- as.data.frame(reMat) %>%
    select(all_of(sel)) %>%
    rename_all( ~ reNames)
  
  re.se <- as.data.frame(reMatSe) %>%
    select(all_of(sel)) %>%
    rename_all( ~ reNames)
  
  
  if (to == "table") {
    
    if (ifSE) {
      reTab <- mapply(function(mean, se) {
        mean <- formatC(mean, format = "f", digits = 3)
        se   <- formatC(se, format = "f", digits = 3)
        paste0(mean, "(", se, ")")
      }, mean = re, se = re.se)
      rownames(reTab) <- rownames(re.se)
    } else {
      reTab <- re
    }
    
    reTab <- tibble::rownames_to_column(as.data.frame(reTab), var = "Data")
    openxlsx::write.xlsx(reTab, file = paste0("../results/ss_", target, ".xlsx"), colNames = T)
    print(reTab)
  }
  
  if (to == "figure") {
    
    # to figures
    plotDat <- as.data.frame(re) %>% 
      rownames_to_column(var = "Data")  %>%
      pivot_longer(-1, values_to = "Mean") %>%
      left_join(as.data.frame(re.se) %>% 
                  rownames_to_column(var = "Data")  %>%
                  pivot_longer(-1, values_to = "SE"), by = c("Data", "name")) %>%
      mutate(method = ifelse(grepl("Stk\\(SS", name), "Stk(SS)", "Stk(PF)"),
             name = factor(name, levels = colnames(re), labels = colnames(re))) %>%
      # add filter
      filter(grepl(",NLL", name))
    
    plotList <- list()
    for(i in unique(plotDat$Data)) {
      plotList[[i]] <- ggplot(plotDat[plotDat$Data == i,], aes(x = name, y = Mean, color = method))+ 
        geom_line(group = 1, color = 1) + 
        geom_pointrange(aes(ymin = Mean - SE, ymax = Mean + SE)) +
        # geom_text(aes(label = round(Mean, 3)), vjust = -0.5) +
        theme_bw() + 
        labs(title = i,
             x = "",
             y = case_when(target == "aucs" ~ "iAUC",
                           target == "ibss" ~ "IBS",
                           target == "logLiks"  ~ "iNNL")) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.title = element_blank(),
              legend.position = "none")
    }
    
    ggpubr::ggarrange(plotlist = plotList, nrow = 2, ncol = 5, common.legend = FALSE)
  }
}

vsSS(target = "aucs", ifSE = TRUE, to = "table") 
vsSS(target = "ibss", ifSE = TRUE, to = "table") 
vsSS(target = "logLiks", ifSE = TRUE, to = "table")

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Show the distribution of weights
#-------------------------------------------------------------------------------

# IPFLasso, PF distribution
pfs <- list()

for(nam in datNames) {
  pf <- c()
  for(i in 1:(nfolds * nsim)){
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    load(paste0("../output/basic_", nam, "_", samID, "_", j, ".RData"))
    pf <- c(pf, paste0(mod$ipf.las$pf, collapse = "-"))
  }
  
  tmp <- setNames(rep(0, 4), c("1-1", "1-2", "1-4", "1-8"))
  pf  <- table(pf) / length(pf)
  tmp[names(pf)] <- pf
  pfs[[nam]] <- tmp
}

# data processing
pfs <- as.data.frame(pfs) %>%
  rownames_to_column(var = "PF") %>%
  pivot_longer(cols = -1,
               names_to = "Data", 
               values_to = "PF(%)") %>%
  mutate(PF = paste0("(", gsub("-", ":", PF), ")"))

gBar <- ggplot(pfs, aes(x = factor(PF), y = `PF(%)`)) + 
  geom_bar(stat="identity") +
  theme_bw() + 
  labs(x = NULL, title = "PF in IPFLasso") +
  facet_grid(~ Data) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 

# functions to get weights of stacking
getAlpha <- function(datName, whichMethod, whichLoss, paraList) {
  
  nfolds = paraList$nfolds
  nsim   = paraList$nsim
  
  alpha <- vector("list", length(whichLoss))
  
  for(i in 1:(nfolds * nsim)){
    
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    
    # merge data
    # get stk
    whichRun = "stk"
    load(paste0("../output/stk_", datName, "_", samID, "_", j, ".RData"))
    
    mod0 <- mod
    # get stk.subset 
    load(paste0("../output/stkPF.pro_", datName, "_", samID, "_", j, ".RData"))
    mod <- c(mod0, mod)
    rm(mod0)
    
    a <- mod[[whichMethod]]$alpha[whichLoss]
    
    for(k in seq_along(a)) {
      if (names(a)[k] %in% c("logLik:lasso", "logLik:ridge"))
        sub_a <- a[[k]][-1]
      else
        sub_a <- a[[k]]
      
      alpha[[k]] <- rbind(alpha[[k]], sub_a)
    }
  }
  names(alpha) = whichLoss
  return(alpha)
}

alphas <- list()
for(datName in datNames){ 
  print(datName)
  alphas[[datName]] <- getAlpha(datName, 
                                whichMethod = "stk.pf", 
                                whichLoss = c("ibsLoss", "logLik:ridge"),
                                paraList = paraList)
}

# Make boxplot 
# derive the data used in ggplot

toggPlotData <- function(alpha, alphaName) {
  re <- c()
  for(i in seq_along(alpha)){
    a <- alpha[[i]]
    a <- a %>% 
      as.data.frame() %>%
      setNames(alphaName) %>%
      pivot_longer(cols = everything(.), 
                   names_to = "subModels",
                   values_to  = "w") %>%
      mutate(lossType = names(alpha)[i])
    
    re <- rbind(re, a)
  }
  return(re)
}

# get data
ggData <- c()
for(nam in names(alphas)) {
  alpha <- alphas[[nam]]
  tmp <- toggPlotData(alpha, c("PF(1:1)", "PF(1:2)", "PF(1:4)", "PF(1:8)"))
  tmp$data <- nam
  ggData <- rbind(ggData, tmp)
}

# boxplot and barPlot
gBox1 <- ggplot(ggData[ggData$lossType == "ibsLoss",], aes(x = factor(subModels), y = w)) + 
  geom_boxplot(outliers = FALSE) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL, title = "PF weights in Stk(IBS)") + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  facet_grid( ~ data) 
# theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10, unit = "mm"))

gBox2 <- ggplot(ggData[ggData$lossType == "logLik:ridge",], aes(x = factor(subModels), y = w)) + 
  geom_boxplot(outliers = FALSE) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL, title = "PF weights in Stk(NLL)") + 
  facet_grid( ~ data) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

gg <- ggpubr::ggarrange(gBar, gBox1, gBox2, nrow = 3)
# output pdf and jpeg format
ggsave("../results/pf_weights.pdf", plot = gg, width = 10, height = 8)
ggsave("../results/pf_weights.jpeg", plot = gg, width = 10, height = 8, units = "in", dpi = 300)

# This plot was used in my presentation: clinical and omics data

# alphas <- list()
# for(datName in datNames){ 
#   print(datName)
#   alphas[[datName]] <- getAlpha(datName, 
#                                 whichMethod = "stk.comp", 
#                                 whichLoss = c("ibsLoss", "logLik:lasso"),
#                                 paraList = paraList)
# }
# 
# # get data
# ggData <- c()
# for(nam in names(alphas)) {
#   alpha <- alphas[[nam]]
#   tmp <- toggPlotData(alpha, c("Clin", "Mol"))
#   tmp$data <- nam
#   ggData <- rbind(ggData, tmp)
# }
# 
# gBox1 <- ggplot(ggData[ggData$lossType == "ibsLoss",], aes(x = factor(subModels), y = w)) + 
#   geom_boxplot(outliers = FALSE) + 
#   theme_bw() + 
#   labs(y = "Weights", x = NULL, title = "CP weights in Stk(IBS)") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
#   facet_grid( ~ data) 
# 
# ggsave("../results/cp_alpha.jpeg", plot = gBox1, width = 8, height = 2.5, units = "in", dpi = 300)

################################################################################
# The distribution of weights in Subset Stacking
################################################################################

alphas <- list()
for(datName in datNames) {
  alphas[[datName]] <- getAlpha(datName, 
                                whichMethod = "stk.subset", 
                                whichLoss = c("ibsLoss", "logLik:ridge"),
                                paraList = paraList)
}

# Boxplot 
toggPlotData <- function(alpha, alphaName) {
  
  re <- c()
  for(i in seq_along(alpha)){
    a <- alpha[[i]]
    a <- a %>% 
      as.data.frame() %>%
      setNames(alphaName) %>%
      pivot_longer(cols = everything(.), 
                   names_to = "subModels",
                   values_to  = "w") %>%
      mutate(lossType = names(alpha)[i])
    
    re <- rbind(re, a)
  }
  return(re)
}

# get data for ggplot
ggData <- c()
for(nam in names(alphas)) {
  alpha <- alphas[[nam]]
  subModName <- paste0("S", rep(1:10, each = 4), ",",
                       rep(c("1-1", "1-2", "1-4", "1-8"), times = 10)) 
  tmp <- toggPlotData(alpha, subModName) %>%
    mutate(data = nam,
           PF = sub(".*\\,(.*)", "\\1", subModels),
           SS = sub("\\,.*$", "", subModels),
           SS = factor(SS, levels = paste0("S", 1:10)),
           subModels = factor(subModels, levels = subModName)) 
  ggData <- rbind(ggData, tmp)
}

# boxplot
# the distribution of weights based on ibs loss
gBox1 <- ggplot(ggData[ggData$lossType == "ibsLoss", ], 
                aes(x = factor(subModels), y = w, fill = SS)) + 
  geom_boxplot(outliers = TRUE,  outlier.size = 0.5) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL) + 
  theme(axis.text.x = element_text(
    angle = 90, vjust = 0.5, hjust = 1, size = 6),
    legend.title = element_blank()) + 
  facet_wrap(vars(data), nrow = 5)

ggsave("../results/ss_weight_ibs.jpeg", gBox1, width = 8, height = 8, units = "in", dpi = 300)

# the distribution of wieghts based on NLL loss
gBox2 <- ggplot(ggData[ggData$lossType == "logLik:ridge",], 
                aes(x = factor(subModels), y = w, fill = SS)) + 
  geom_boxplot(outliers = TRUE, outlier.size = 0.5) + 
  theme_bw() + 
  labs(y = "Weights", x = NULL) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                   size = 6),
        legend.title = element_blank()) + 
  facet_wrap(vars(data), nrow = 5) 

ggsave("../results/ss_weight_logLik.jpeg", gBox2, width = 8, height = 8, units = "in", dpi = 300)

#-------------------------------------------------------------------------------
# Computational complexity
#-------------------------------------------------------------------------------

# get time Data

getTime <- function(datName, paraList){
  
  nfolds = paraList$nfolds
  nsim = paraList$nsim
  
  # select the method
  methodNames <- c("Clin", "Mol", "Naive", "FS(las)", "PCR(las)", "IPFLasso", "Stk(PF)", "Stk(CP)", "Stk(SS)")
  
  time <- c()
  for(i in 1:(nfolds * nsim)){
    samID <- ceiling(i / nfolds)
    j <- i - (samID - 1) * nfolds
    
    mods <- list()
    load(paste0("../output/basic_", datName, "_", samID, "_", j, ".RData"))
    mods <- c(mods, mod)
    load(paste0("../output/stk_", datName, "_", samID, "_", j, ".RData"))
    mods <- c(mods, mod[-1])
    load(paste0("../output/stkPF.pro_", datName, "_", samID, "_", j, ".RData"))
    mods <- c(mods, mod[-1])
    
    # mods <- mods[-length(mods)]
    mod <- mods 
    rm(mods)
    
    tmp <- lapply(mod, `[[`, "time")
    time <- rbind(time, unlist(lapply(tmp, as.numeric, units = "mins"))) # be careful of this one
  }
  
  # get data for ggplot
  timeggData <- time %>% 
    as.data.frame() %>%
    setNames(methodNames) %>%
    pivot_longer(everything(),
                 names_to = "Method",
                 values_to = "Time") %>%
    group_by(Method) %>%
    summarise(Time = mean(Time)) %>% 
    mutate(Method = factor(Method, levels = methodNames))
  
  return(timeggData)
}


#- using BLCA as an example, ggplot

timeggData <- getTime("BLCA", paraList) %>%  
  arrange(Time) %>%
  mutate(Method = factor(Method, levels = Method))

p1 <- ggplot(timeggData, aes(x = Method, y = Time)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 200)) +
  labs(x = NULL, title = "A. BLCA", y = "Time, minutes") + 
  coord_flip()

timeggData <- getTime("BRCA", paraList)  %>%  
  arrange(Time) %>% 
  mutate(Method = factor(Method, levels = Method))

p2 <- ggplot(timeggData, aes(x = Method, y = Time)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_bw() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 700)) +
  labs(x = NULL, title = "B. BRCA", y = "Time, minutes") + 
  coord_flip()

gg <- ggpubr::ggarrange(p1, p2, nrow = 1)

# 
# output figures in pdf and jpeg format
ggsave("../results/timeCost.pdf", plot = gg, width = 7, height = 3)
ggsave("../results/timeCost.jpeg", plot = gg, width = 7, height = 3, units = "in", dpi = 600)


# stop here

#-------------------------------------------------------------------------------






