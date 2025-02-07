
# load the package 

loadPackages <- function(package_list) {
    for (package in package_list) {
        if (require(package, character.only = TRUE, quietly = TRUE)) {
            # Package is installed, load it
            message(paste("Package", package, "loaded successfully."))
        } else {
            # Package is NOT installed, install and then load
            message(paste("Package", package, "not found. Installing..."))
            tryCatch({
                install.packages(package, dependencies = TRUE)  # Install with dependencies
                if (require(package, character.only = TRUE)) {
                    message(paste("Package", package, "installed and loaded successfully."))
                } else {
                    stop(paste("Failed to load package", package, "after installation."))
                }
            }, error = function(e) {
                message(paste("Error installing package", package, ":", e$message))
                # You could choose to stop the script here with stop() if a critical package fails.
            })
        }
    }
    invisible(NULL)
}

#- function to import data

getData <- function(datName, log2 = T, toBinary = T, exclude = ""){
    
    # datName, data name
    # log2, whether log2 transformation on TPM (molecular information)
    # exclude, exclude a variable in X
    
    datDir  <- paste0("../data/", datName, ".RData")
    load(datDir)
    
    y <- survival::Surv(X[,"time"], event = X[,"status"])
    X <- X[,-c(1:3)]
    
    if (exclude != "") {
        
        if (!all(exclude %in% colnames(X))){
            warning(paste0("exclude variables not in the data:", 
                           paste0(setdiff(exclude, colnames(X)), collapse = ",", ".")))
        }
        
        ind <- na.omit(match(exclude, colnames(X)))
        if (length(ind) > 0)
            X <- X[ , -ind]
    }
    
    b <- as.numeric(!grepl("clin", colnames(X))) + 1
    
    if (log2) {
        X[,b == 2] <- log2(X[,b == 2] + 0.001)
    }
    
    #- trans to binary response
    if (toBinary){
        y <- biY(y)
        X <- X[y$sel,]
        y <- y$biY
    }
    
    #- Rename molecules since the original name is too long
    colnames(X) <- c(colnames(X)[b == 1], 
                     paste0("m", seq_len(sum(b == 2))))
    
    #- return
    .GlobalEnv$X <- as.matrix(X)
    .GlobalEnv$y <- y
    .GlobalEnv$blocks <- b
}



getSplits <- function(n, nsim, nfolds, seed, foldid.internal = T){
    
    # n, sample size
    # nsim, number of duplicated runs
    
    set.seed(seed)
    
    sampleSplits <- list()
    foldids <- list()
    
    for(i in seq_len(nsim)){
        samsplit <- sample(rep(1:10, length = n), n)
        sampleSplits[[i]] <- samsplit
        foldids[[i]] <- sapply(1:nfolds, 
                               function(id){ 
                                   nlength = sum(samsplit != id)
                                   sample(rep(1:10, length = nlength), nlength)}, 
                               simplify = F)
    }
    
    if (nsim == 1) {
        sampleSplits <- unlist(sampleSplits)
        foldids <- foldids[[1]]
    }
    
    if (foldid.internal) {
        
        foldid.internals <- list()
        for(i in seq_len(nsim)){
            tmp <- list()
            for(j in seq_along(foldids[[i]])){
                sams <- foldids[[i]][[j]]
                tmp[[j]] <- sapply(1:nfolds, 
                                   function(x){ 
                                       nlength = sum(sams != x)
                                       sample(rep(1:10, length = nlength), nlength)
                                   }, simplify = F)
            }
            foldid.internals[[i]] <- tmp
        }
        .GlobalEnv$foldid.internals <- foldid.internals
    }
    
    #- return
    .GlobalEnv$sampleSplits <- sampleSplits
    .GlobalEnv$foldids <- foldids
}

#-

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

#-------------------------------------------------------------------------------
#- prediction performance with deviance difference
#-------------------------------------------------------------------------------
#- predict the survival P of given times

predProbs <- function(X.v, obj, times, ifCVMods = FALSE){
    
    #- clin
    getLp <- function(mod, X.v, times){
        
        #- mod, subset of obj
        beta <- mod$beta
        if (inherits(mod$beta, "dgCMatrix"))
            s <- seq_along(beta)
        else
            s <- names(beta)
        
        # all coefficients are 0, if no variable selected
        if (length(s) == 0)  {
            s = seq_len(ncol(X.v))
            beta = rep(0, ncol(X.v))
        }
        
        lp.v <- X.v[,s, drop = F] %*% beta
        chz <- cumhz(mod$lp, mod$y, times)
        
        # metaSurv <- matrix(NA, nrow = nrow(lp.v), ncol = length(times),
        #                 dimnames = list(rownames(X.v), 
        #                                 as.character(times)))
        
        metaSurv <- c()
        for(j in seq_along(times)){
            metaSurv <-  cbind(metaSurv, drop(exp(-chz[j] * exp(lp.v))))
        }
        colnames(metaSurv) <- as.character(times)
        
        metaSurv
    }
    
    ind <- !grepl("stk", names(obj))
    re <- mapply(getLp, mod = obj[ind], MoreArgs = list(X.v = X.v, times = times), 
                 SIMPLIFY = F)
    
    #- stacking -pf
    
    stkMod <- names(obj)[grepl("stk", names(obj))]
    
    for(i in stkMod) {
        
        if (grepl("super", i)) {
            # #- stacking -combine clin, mol, naive, pca, ipf
            re[[i]] <- predict.stack.super.cox(obj[[i]], X.v, times = times)
        } else {
            re[[i]] <- predict.stack.cox(obj[[i]], X.v, times = times, ifCVMods = ifCVMods)
        }
    }
    re
}

#-------------------------------------------------------------------------------
# Predicton Measure
#-------------------------------------------------------------------------------

# to calculate the IBS with given time
getTimeIBS <- function(survProb, 
                       Xclin.t, Xclin.v, y.t, y.v, 
                       times, cens.method,
                       sumRes = TRUE,
                       truncICPW = 1){
    
    # matrix, predicted survival probabilities of the validation data
    # y.t, response of training data
    # y.v, response of validation data
    # times to integrate the BS, corresponding to the columns of survProb
    # cens.method, the methods to calcualte censoring weights. "Breslow", "KM", "equal".
    # sumRes, summarize the result 
    
    # ibsDat <- matrix(NA, length(times), 2, dimnames = list(seq_along(times), c("time", "ibs")))
    
    bss <- c()
    for(i in seq_along(times)){
        bs <- brierScore(Xclin.t, Xclin.v, y.t, y.v, survProb[,i], times[i], cens.method = cens.method)
        bs[,"time"] <- times[i]
        bss <- rbind(bss, bs)
    }
    
    # truncation on weights
    bss[,1] <- ifelse(bss[,1] >= quantile(bss[,1], truncICPW), quantile(bss[,1], truncICPW), bss[,1])
    
    ibsDat <- bss %>% 
        group_by(time) %>%
        summarise(ibs = sum(weight*(Zind - survProb)^2) / sum(weight)) %>%
        as.matrix()
        
    # truncation
    # for(i in seq_along(times)){
    #     bs <- brierScore(Xclin.t, Xclin.v, y.t, y.v, survProb[,i], times[i], cens.method = cens.method)
    #     ibsDat[i,"time"] <- times[i]
    #     ibsDat[i,"ibs"]  <- sum(bs[,"weight"] * (bs[,"Zind"] - bs[,"survProb"]) ^ 2) / sum(bs[,"weight"])
    # }
    # 
    if (sumRes) {
        ibsDat <- na.omit(ibsDat)
        tRange <- abs(range(ibsDat[,"time"])[1] - range(ibsDat[,"time"])[2])
        trapezoidal(ibsDat[,"time"], ibsDat[,"ibs"]) / tRange
    } else {
        ibsDat
    }
    
}


# to calculate AUC with given time
getTimeAUC <- function(survProb, Xclin.t, Xclin.v, y.t, y.v, 
                       times, cens.method = "equal", truncICPW = 1, 
                       sumRes = TRUE) {
    
    # matrix, predicted survival probabilities of the validation data
    # y.t, response of training data
    # y.v, response of validation data
    # times to integrate the BS, corresponding to the columns of survProb
    # cens.method, the methods to calcualte censoring weights. "Breslow", "KM", "equal".
    # derive from timeROC function in timeROC package
    # https://github.com/cran/timeROC/blob/master/R/timeROC_3.R#L106 
    
    calAUC <- function(marker, events, weight){
        
        ordMarker <- order(marker, decreasing = T)
        events    <- events[ordMarker]
        weights   <- weight[ordMarker]
        marker    <- marker[ordMarker]
        
        # TP
        TP   <- c(0, cumsum(weights * events)) /  sum(weights * events)
        TP   <- TP[!duplicated(marker[ordMarker])]
        
        # FP
        FP   <- c(0, cumsum(weights * (1 - events))) /  sum(weights * (1 - events))
        FP   <- FP[!duplicated(marker[ordMarker])]
        
        AUC <- trapezoidal(FP, TP)
        AUC
    }
    
    #-
    aucs <- matrix(NA, nrow = length(times), ncol = 2, dimnames = list(seq_along(times), c("time", "auc")))
    # bss[,1] <- ifelse(bss[,1] >= quantile(bss[,1], truncICPW), quantile(bss[,1], truncICPW), bss[,1])
    
    for(i in seq_along(times)){
        
        bs <- brierScore(Xclin.t, Xclin.v, y.t, y.v, survProb[,i], times[i], cens.method = cens.method)
        bs[,1] <- ifelse(bs[,1] >= quantile(bs[,1], truncICPW), quantile(bs[,1], truncICPW), bs[,1])
        aucs[i,"auc"]  <- calAUC(marker = bs[,3], events = bs[,2], weight = bs[,1])
        aucs[i,"time"] <- times[i]
    }
    
    if (!sumRes) {
        aucs
    } else {
        aucs <- na.omit(aucs)
        tRange <- abs(range(aucs[,"time"])[1] - range(aucs[,"time"])[2])
        trapezoidal(aucs[,"time"], aucs[,"auc"]) / tRange
    }
}


getTimelogLik <- function(survProb, y.v, times, sumRes = TRUE) {
    
    time_i  <- y.v[,1]
    event_i <- y.v[,2]
    logDat <- c()
    
    for(i in seq_along(times)){
        
        tmpZind <- ifelse(time_i > times[i], 1, 0 / event_i)
        tmpDat  <- data.frame(survP = survProb[ ,i],
                              timeObs = time_i, 
                              event = event_i, 
                              Zind = tmpZind, 
                              time = times[i])
        logDat <- rbind(logDat, tmpDat)
    }
    
    rmlogDat <- logDat[ ,"Zind"] %in% c(NaN, Inf, NA)
    logDat <- logDat[!rmlogDat, ] %>%
        as.data.frame() %>%
        group_by(time) %>%
        summarise(logLik = sum(-1 * {Zind * log(survP) + (1 - Zind) * log(1 - survP)})) %>%
        as.matrix()
    
    if (!sumRes) {
        logDat
    } else {
        logDat <- na.omit(logDat)
        tRange <- abs(range(logDat[,"time"])[1] - range(logDat[,"time"])[2])
        trapezoidal(logDat[,"time"], logDat[,"logLik"]) / tRange
        # sum(logDat[,2])
    }
}

# time dependent AUC of ROC using SurvProb
# trapezoidal rule

trapezoidal <- function(x, y){
    
    ord  <- order(x, decreasing = F)
    y   <- y[ord]
    x   <- x[ord]
    
    n    <- length(y)
    disty <- (y[-1] + y[-n]) / 2
    distx <- x[-1] - x[-n]
    
    sum(disty * distx)
}

#-------------------------------------------------------------------------------

