
#-------------------------------------------------------------------------------
#- functions to fit models 
#- fit.basic: clinical model, molecular model, naive model, dimensional reduction based on PCA, dimension reduction based on fs, IPFLasso
#- fit.stk includes stk(PF) and stk(CP)
#- fit.stkPF.pro includes subset stacking
#-------------------------------------------------------------------------------

fit.basic <- function(X, y, 
                      ncomp = NULL,
                      blocks = rep(1, ncol(X)),
                      foldid,
                      times,
                      cox.time,
                      intercept = T, 
                      foldid.internal = NULL, ...) {
  
  require(glmnet)
  
  # clinical model
  t1 <- Sys.time()
  s  <- blocks == 1
  cvs  <- cv.glmnet(X[,s], y, family = "cox", foldid = foldid)
  clin <- glmnet(X[,s], y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()
  
  beta <- coef(clin)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  clin <- list(beta = beta,
               lp   = lp,
               y    = y,
               time = t2 - t1)
  
  #- molecular model
  t1 <- Sys.time()
  s  <- blocks == 2
  cvs <- cv.glmnet(X[,s], y, family = "cox", foldid = foldid)
  mol <- glmnet(X[,s], y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()
  
  beta <- coef(mol)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  mol <- list(beta = beta,
              lp   = lp,
              y    = y,
              time = t2 - t1)
  
  #- naive model
  t1 <- Sys.time()
  cvs <- cv.glmnet(X, y, family = "cox", foldid = foldid)
  naive <- glmnet(X, y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()
  
  lp   <- predict(naive, X)
  beta <- coef(naive)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  naive <- list(beta = beta,
                lp   = lp,
                y    =y,
                time = t2 - t1)
  
  #- dimension reduction based on feature screening
  t1 <- Sys.time()
  cv.fs <- cv.drfs(X, y,
                   family = "cox",
                   foldid = foldid,
                   nfolds = 10,
                   blocks = blocks,
                   blocks.screen = 2,
                   penalty.factors = list(c(1, 1)),
                   cuts = seq(0.1, 0.5, 0.1),
                   alpha = 1)
  
  fs.las <- drfs(X, y,
                 family = "cox",
                 blocks = blocks,
                 blocks.screen = 2,
                 feature.order = cv.fs$feature.order,
                 cut = cv.fs$cut.min,
                 pf  = cv.fs$pf.min,
                 alpha  = 1,
                 lambda = cv.fs$lambda.min)
  
  t2 <- Sys.time()
  
  beta <- coef(fs.las)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  fs.las <- list(beta = beta,
                 lp   = lp,
                 y    = y,
                 cut = cv.fs$cut.min,
                 feature.cut = fs.las$feature.cut,
                 time = t2 - t1)
  
  # dimension reduction based on PCA
  if (is.null(ncomp))
    ncomp = nrow(X)
  
  t1 <- Sys.time()
  cv.pcr <- cv.drpcr(X, y,
                     alpha = 1,
                     family = "cox",
                     ncomp = ncomp,
                     blocks = blocks,
                     blocks.pc = 2,
                     penalty.factors = list(c(1, 1)),
                     foldid = foldid)
  
  pcr.las <- drpcr(X, y,
                   alpha = 1,
                   ncomp = ncomp,
                   family = "cox",
                   blocks = blocks,
                   blocks.pc = 2,
                   pf = cv.pcr$pf.min,
                   lambda = cv.pcr$lambda.min)
  t2 <- Sys.time()
  
  beta = pcr.las$back.mod$beta
  lp   = getLp(X, beta)
  
  pcr.las <- list(beta = beta,
                  lp   = lp,
                  y    = y,
                  beta.pcr = c(pcr.las$a0, pcr.las$beta[,1]),
                  centers = pcr.las$centers,
                  scales  = pcr.las$scales,
                  blocks = pcr.las$blocks,
                  time =  t2 - t1)
  
  #- ipflasso
  t1 <- Sys.time()
  cv.ipf <- cv.ipflas(X, y,
                      family = "cox",
                      foldid = foldid,
                      blocks = blocks,
                      alpha = 1,
                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)))
  
  ipf.las <- ipflas(X, y,
                    family = "cox",
                    alpha = 1,
                    blocks = blocks,
                    pf = cv.ipf$pf.min,
                    lambda = cv.ipf$lambda.min)
  t2 <- Sys.time()
  
  beta <- coef(ipf.las)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  ipf.las <- list(beta = beta,
                  lp   = lp,
                  y    = y,
                  pf   = cv.ipf$pf.min,
                  time = t2 - t1)
  
  # return a list
  list(clin = clin,
       mol = mol,
       naive = naive,
       fs.las = fs.las,
       pcr.las = pcr.las,
       ipf.las = ipf.las)
}


# stacking ---------------------------------------------------------------------

fit.stk <- function(X, y, 
                    blocks = rep(1, ncol(X)),
                    foldid,
                    times,
                    cox.time,
                    intercept = TRUE, 
                    foldid.internal = NULL, ...) {
  
  require(glmnet)
  
  #- clinical models, for validation
  t1 <- Sys.time()
  s  <- blocks == 1
  cvs  <- cv.glmnet(X[,s], y, family = "cox", foldid = foldid)
  clin <- glmnet(X[,s], y, family = "cox", lambda = cvs$lambda.min)
  t2 <- Sys.time()
  
  beta <- coef(clin)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  clin <- list(beta = beta,
               lp   = lp,
               y    = y,
               time = t2 - t1)
  
  # stacking based on penalty factor
  t1 <- Sys.time()
  stk.pf <- stack.cox(X, y,
                      times = times,
                      foldid = foldid,
                      foldid.internal = foldid.internal,
                      nfold = 10,
                      blocks = blocks,
                      block.clin = 1,
                      ifCox = TRUE, 
                      cox.time = cox.time,
                      method = "bypf",
                      cens.method = "clinicalAdjusted",
                      penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                      track = TRUE, 
                      optimLoss = c("ibsLoss", "logLik", "PaLogLik"), 
                      optimFun = c("lasso", "ridge"),
                      moreOptimFun = NULL,
                      intercept = TRUE,
                      foldidSL = NULL,
                      foldidSL.cox = NULL,
                      nfoldSL = 10,
                      limits = 0,
                      ifWeights = c(TRUE, TRUE, FALSE),
                      truncICPW = 0.95, 
                      truncValue = 0.01, 
                      transPFun = function(x) log(x / (1 - x)))
  t2 <- Sys.time()
  stk.pf$time = t2 - t1
  
  #- stacking based on component-wise sub-models
  t1 <- Sys.time()
  stk.comp <- stack.cox(X, y,
                        times = times,
                        foldid = foldid,
                        foldid.internal = foldid.internal,
                        nfold = 10,
                        blocks = blocks,
                        block.clin = 1,
                        ifCox = TRUE, 
                        cox.time = cox.time,
                        method = "bycomp",
                        cens.method = "clinicalAdjusted",
                        penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                        track = TRUE, 
                        optimLoss = c("ibsLoss", "logLik", "PaLogLik"), 
                        optimFun = c("lasso", "ridge"),
                        moreOptimFun = NULL,
                        intercept = TRUE,
                        foldidSL = NULL,
                        foldidSL.cox = NULL,
                        nfoldSL = 10,
                        limits = 0,
                        ifWeights = c(TRUE, TRUE, FALSE),
                        truncICPW = 0.95, 
                        truncValue = 0.01, 
                        transPFun = function(x) log(x / (1 - x)))
  
  t2 <- Sys.time()
  stk.comp$time = t2 - t1
  
  list(clin = clin, 
       stk.pf = stk.pf,
       stk.comp = stk.comp)
}

# subset stacking, extension of stk(pf)
fit.stkPF.pro <- function(X, y, 
                          blocks = rep(1, ncol(X)),
                          foldid,
                          times,
                          cox.time,
                          intercept = TRUE, 
                          foldid.internal = NULL, ...) {
  
  require(glmnet)
  
  #- clinical models, for validation
  t1 <- Sys.time()
  s  <- blocks == 1
  cvs  <- cv.glmnet(X[,s], y, family = "cox", foldid = foldid)
  clin <- glmnet(X[,s], y, family = "cox", lambda = cvs$lambda.min)
  
  beta <- coef(clin)[,1]
  beta <- beta[beta != 0]
  lp   <- getLp(X, beta)
  
  t2 <- Sys.time()
  
  clin <- list(beta = beta,
               lp   = lp,
               y    = y,
               time = t2 - t1)
  
  
  #- subset stacking
  t1 <- Sys.time()
  indClin <- blocks == 1
  blocks[!indClin] <- sort(rep(2:11, length = sum(!indClin)))
  
  for(k in 2:11) {
    
    ind <- blocks %in% c(1, k)
    b   <- as.numeric(blocks[ind] != 1) + 1
    stkTmpt <- StkCoxBL(X[ ,ind], y, 
                        times = times,
                        foldid = foldid,
                        foldid.internal = foldid.internal,
                        blocks = b, 
                        block.clin = 1,
                        nfold = 10,
                        ifCox = TRUE, 
                        cox.time = cox.time,
                        method = "bypf",
                        cens.method = "clinicalAdjusted",
                        penalty.factors = list(c(1, 1), c(1, 2), c(1, 4), c(1, 8)),
                        track = TRUE)
    
    if (k == 2) {
      stk.subset <- stkTmpt
    } else {
      stk.subset <- merge.stk_cox(stk.subset, stkTmpt)
    }
  }
  # 
  subNames <- paste0("sub", rep(1:10, each = 4), "_", c("1-1", "1-2", "1-4", "1-8"))
  
  names(stk.subset$subMod$beta)  <- subNames
  names(stk.subset$subMod$lp)    <- subNames
  colnames(stk.subset$CVpredCox) <- subNames
  colnames(stk.subset$CVpred) <- c("weight", "Zind", subNames)


  for(j in seq_along(stk.subset$CVMods)) {
    names(stk.subset$CVMods[[j]]$mods) <- subNames
  }
  
  # super learner
  stk.subset <- StkCoxSL(stk.subset, 
                         optimLoss = c("ibsLoss", "logLik", "PaLogLik"), 
                         optimFun = c("ridge"), 
                         moreOptimFun = NULL, 
                         intercept = TRUE, 
                         foldid = NULL, foldid.cox = NULL, nfold = 10, 
                         limits = 0, 
                         ifWeights = c(TRUE, TRUE, FALSE), 
                         truncICPW = 0.95, 
                         truncValue = 0.01, 
                         transPFun = function(x) log(x / (1 - x))) 
  
  t2 <- Sys.time()
  stk.subset$time <- t2 - t1
  
  return(list(clin = clin, stk.subset = stk.subset))
}

#-------------------------------------------------------------------------------



















