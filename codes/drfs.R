

coxph.fit.p <- function(x, y, strata, offset, init, control,
                        weights, method, nocenter = NULL) {
  
  n <- nrow(y)
  if (is.matrix(x)) {
    nvar <- ncol(x) 
  } else {
    if (length(x)==0) {
      nvar <-0
    } else {
      nvar <- 1
      x <- matrix(x, ncol = 1)
    }
  }
  time <- y[,1]
  status <- y[,2]
  
  # Sort the data (or rather, get a list of sorted indices)
  if (missing(strata)) {
    sorted <- order(time)
    strata <- NULL
    newstrat <- as.integer(rep(0,n))
  } else {
    sorted <- order(strata, time)
    strata <- strata[sorted]
    newstrat <- as.integer(c(1*(diff(as.numeric(strata))!=0), 1))
  }
  
  if (missing(offset) || is.null(offset)) offset <- rep(0,n)
  if (missing(weights)|| is.null(weights)) {
    weights<- rep(1, n)
  } else {
    if (any(weights<=0)) stop("Invalid weights, must be >0")
    weights <- weights[sorted]
  }
  stime <- as.double(time[sorted])
  sstat <- as.integer(status[sorted])
  
  if (nvar==0) {
    # A special case: Null model.
    #  (This is why I need the rownames arg- can't use x' names)
    # Set things up for 0 iterations on a dummy variable
    x <- as.matrix(rep(1.0, n))
    nullmodel <- TRUE
    nvar <- 1
    init <- 0
    maxiter <- 0
  } else {
    nullmodel <- FALSE
    maxiter <- control$iter.max
    if (!missing(init) && length(init)>0) {
      if (length(init) != nvar) stop("Wrong length for inital values")
    }
    else init <- rep(0, nvar)
  }
  
  # 2012 change: individually choose which variable to rescale
  # default: leave 0/1 variables along
  if (is.null(nocenter)) zero.one <- rep(FALSE, ncol(x))
  else zero.one <- apply(x, 2, function(z) all(z %in% nocenter))
  
  storage.mode(weights) <- storage.mode(init) <- "double"
  coxfit <- .Call(survival:::Ccoxfit6, 
                  as.integer(maxiter),
                  stime, 
                  sstat,
                  x[sorted,],
                  as.double(offset[sorted]),
                  weights,
                  newstrat,
                  as.integer(method == method),
                  as.double(control$eps),
                  as.double(control$toler.chol),
                  as.vector(init),
                  ifelse(zero.one, 0L, 1L))
  
  beta <- coxfit$coef
  df <- length(beta[!(is.na(beta))])
  logtest <- -2 * (coxfit$loglik[1] - coxfit$loglik[2])
  pvalue = pchisq(logtest, df, lower.tail = FALSE)

  return(pvalue)
}          

#- dimension reduction with feature screening ----------------------------------
#' fs 

#' @description 
#' Feature screening method using coxPH for survival and glm for logistic regression
#' 
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' 
#' @examples 
#' 
#' @import survival
#' @export fs


fs <- function(X, y, family, offset = NULL, weights = rep(1, nrow(X))){
  
  n <- length(y)
  
  if (is.null(offset))
    offset = rep.int(0, n)
  
  # update, using quasibinomial() family
  uni.glm <- function(x, y, offset, weights){
    
    obj <- glm.fit(cbind(1, x), y, 
                   family = quasibinomial(), 
                   offset = offset)
    
    dipersion <- sum((obj$weights * obj$residuals^2)[obj$weights > 0])/ obj$df.residual
    
    Qr  <- obj$qr
    coef.p <- obj$coefficients
    var.cf <- sqrt(diag(chol2inv(Qr$qr))) * dipersion
    
    t.values <- coef.p / var.cf
    2 * pnorm(-abs(t.values))[2]
  }
  
  if (family == "cox")
    p <- apply(X, 2, coxph.fit.p, y = y, weights = weights,
               control = coxph.control(), method = "efron", offset = offset)
      
  if (family == "binomial")
    p <- apply(X, 2, uni.glm, y = y, offset = offset, weights = weights)
  p
  
}

#' cv.drfs 
#' 
#' @description 
#' 
#' Cross-validation function for dimension reduction based on feature screening, a method to combine clinical 
#' molecular variables for prediction. Step1, feature screening. Step2, model fitting
#' 
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' @param blocks block of variables
#' @param blocks.screen block variables for feature screening
#' @param cuts the % of top variables included
#' @param foldid list or vector depends on ncv
#' @param nfolds defaulted as 10
#' @param alpha defaulted as 1 
#' @param blocks vector of the same length of ncol(X)
#' @param penalty.factors list of penalty factors 
#' @param type.measure "default" 
#' @param lambda NULL 
#' @param nlam 100 
#' @param pf penalty factor
#' @param lambda.min.ratio 0.01 
#' @param weights 1 
#' 
#' @examples 
#' 
#' @import glmnet
#' @export cv.drfs

cv.drfs <- function(X, y, 
                    family, 
                    foldid = NULL, 
                    nfolds = 10, 
                    blocks = rep(1, ncol(X)), 
                    blocks.screen = 1,
                    offset = NULL,
                    cuts = seq(0.1, 1, 0.1), 
                    alpha = 1, 
                    lambda = NULL,
                    penalty.factors = NULL, 
                    lambda.min.ratio = 0.001,
                    nlam = 100, 
                    weights = rep(1, nrow(X)),
                    type.measure = "default", ...){
  
  ps <- fs(X, y, family = family, offset = offset, weights = weights)
  ord <- order(ps, decreasing = F)
  
  p  <- ncol(X)
  n  <- nrow(X)
  nblocks <- length(unique(blocks))
  
  ls <- list()
  for(c in cuts){
    #- variables that screening not needed
    in.features <- seq_len(p)[!blocks %in% blocks.screen] 
    
    for(s in blocks.screen){
      
      sel   <- ord[ord %in% which(blocks == s)]
      in.f  <- sel[1:floor(quantile(seq_along(sel), c))]
      in.features <- c(in.features, in.f)
    }
    ls[[paste(c)]] <- in.features
  }
  
  #- foldid
  if (is.null(foldid))
    foldid <- sample(rep(1:nfolds, length = n), n)
  
  if (is.null(penalty.factors)) {
    l <- list()
    for (i in 1:nblocks) {
      l[[i]] <- c(1, 2, 4, 8, 16)
    }
    penalty.factors <- do.call("expand.grid", list(... = l))
    penalty.factors <- penalty.factors[apply(penalty.factors,
                                             1, function(x) any(x == 1)), , drop = F]
    colnames(penalty.factors) <- paste0("B", 1:nblocks)
    penalty.factors <- split(penalty.factors, rep(1:nrow(penalty.factors)),
                             ncol(penalty.factors))
  }
  
  pfs <- list()
  for(i in seq_along(penalty.factors)){
    pfs[[i]] <- as.numeric(as.character(factor(blocks, unique(blocks), 
                                               penalty.factors[[i]])))
  }
  
  #- get lambda:using the first cut
  if(is.null(lambda)){
    if (family == "cox")
      lambda <- get.lambda.cox(X[,ls[[which.min(cuts)]]], y, 
                               alpha = alpha, 
                               pfs = lapply(pfs, `[`, ls[[which.min(cuts)]]), 
                               lambda.min.ratio = lambda.min.ratio, 
                               nlam = nlam,
                               weights = weights)
    else
      lambda = get.lambda(X[,ls[[which.min(cuts)]]], y, 
                          family = family,
                          alpha = alpha,
                          nlam = nlam, 
                          lambda.min.ratio = lambda.min.ratio,
                          weights = weights,
                          pfs = lapply(pfs, `[`, ls[[which.min(cuts)]]))
  }
  
  statlist <- list()
  for(i in seq_along(cuts)){
    for(j in seq_along(penalty.factors)){
      
      #- get blocks --- 
      cv <- cv.glmnet(x = X[,ls[[i]]], y = y,
                      family = family,
                      foldid = foldid,
                      weights = weights,
                      type.measure = type.measure,
                      alpha = alpha,
                      penalty.factor = pfs[[j]][ls[[i]]],
                      lambda = lambda)
      
      statlist[[paste0("cut:", cuts[i])]][[paste0("pf:", paste0(penalty.factors[[j]],collapse = "-"))]] <- cv[1:6]
    }
  }
  
  stat <- unlist(statlist, F)
  cvm  <- as.matrix(Reduce(cbind, lapply(stat, `[[`, "cvm")))
  cvsd <- as.matrix(Reduce(cbind, lapply(stat, `[[`, "cvsd")))
  nzero <- as.matrix(Reduce(cbind, lapply(stat, `[[`, "nzero")))
  
  if (type.measure %in% c("auc", "C"))
    cvm <- -cvm
  
  #- get position
  pos <- which(cvm <= min(cvm), arr.ind = T)
  cuts.ind <- rep(cuts, each = length(penalty.factors))
  pf.ind   <- rep(unlist(lapply(penalty.factors, paste0, collapse = "-")), 
                  times = length(cuts))
  
  #- favoring large cuts and pos
  id.min <- pos[order(pos[,2], pos[,1], decreasing = F)[1], ]
  lambda.min <- lambda[id.min[1]]
  cut.min    <- cuts.ind[id.min[2]]
  pf.min     <- pf.ind[id.min[2]]
  nzero.min  <- nzero[id.min[1], id.min[2]]
  
  #- 1sd rule
  cvm.1se <- cvm[id.min[1], id.min[2]] + cvsd[id.min[1], id.min[2]]
  pos     <- which(cvm <= cvm.1se, arr.ind = T)
  id.1se     <- pos[order(pos[,2], pos[,1], decreasing = F)[1], ]
  lambda.1se <- lambda[id.1se[1]]
  cut.1se    <- cuts.ind[id.1se[2]]
  pf.1se     <- pf.ind[id.1se[2]]
  nzero.1se  <- nzero[id.1se[1], id.1se[2]]
  
  #- summary index
  index = matrix(NA, ncol = 6, nrow = 2,
                 dimnames = list(c("min", "1se"),
                                 c("Lambda", "Cut", "PF", "Measure", "SE", "Nonzero")))
  
  index[,1] <- round(c(lambda.min, lambda.1se), 2)
  index[,2] <- c(paste0(cut.min, collapse = "-"), paste0(cut.1se, collapse = "-"))
  index[,3] <- c(pf.min, pf.1se)
  index[,4] <- round(c(cvm[id.min[1], id.min[2]], cvm[id.1se[1], id.1se[2]]), 2)
  index[,5] <- round(c(cvsd[id.min[1], id.min[2]], cvsd[id.1se[1], id.1se[2]]), 2)
  index[,6] <- c(nzero[id.min[1], id.min[2]], nzero[id.1se[1], id.1se[2]])
  
  fit <- list(statlist = statlist,
              cut = cuts,
              lambda.min = lambda.min,
              lambda.1se = lambda.1se,
              cut.min = cut.min,
              cut.1se = cut.1se,
              pf.min  = as.numeric(unlist(strsplit(pf.min, "-"))),
              pf.1se  = as.numeric(unlist(strsplit(pf.1se, "-"))),
              nzero.min = nzero.min,
              nzero.1se = nzero.1se,
              index = index,
              feature.order = order(ps, decreasing = F))
  fit
}

#' drfs 
#' 
#' @description 
#' 
#' Function for dimension reduction based on feature screening
#'  
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' @param blocks block of variables
#' @param blocks.screen block variables for feature screening
#' @param cut the % of top variables included
#' @param feature.order the feature orders from feature screening
#' @param alpha defaulted as 1 
#' @param blocks vector of the same length of ncol(X)
#' @param penalty.factors list of penalty factors 
#' @param type.measure "default" 
#' @param lambda NULL 
#' @param nlam 100 
#' @param pf penalty factor
#' @param lambda.min.ratio 0.01 
#' @param weights 1 
#' 
#' @examples 
#' 
#' @import glmnet
#' @export cv.drfs

drfs <- function(X, y, family, 
                 blocks = rep(1, ncol(X)), 
                 offset = NULL,
                 blocks.screen = 1,
                 feature.order = NULL,
                 cut = NULL, 
                 alpha = 1, 
                 lambda = NULL,
                 pf = rep(1, ncol(X)), 
                 lambda.min.ratio = 0.001,
                 nlam = 100, 
                 weights = rep(1, nrow(X)), ...) {
  
  if (is.null(feature.order)) {
    
    ps <- fs(X, y, family = family, offset = offset, weights = weights)
    rks <- order(ps, decreasing = F)
  }
  
  p <- ncol(X)
  feature.cut <- seq_len(p)[!blocks %in% blocks.screen] 
  for(s in blocks.screen){
    sel <- feature.order[feature.order %in% which(blocks == s)]
    in.f <- sel[1:floor(quantile(1:p, cut))]
    feature.cut <- c(feature.cut, in.f)
  }
  
  pf <- as.numeric(as.character(factor(blocks, unique(blocks), pf)))
  
  if (is.null(lambda)){
    
    if (family == "cox"){
      lambda <- get.lambda.cox(X[,feature.cut], y, 
                                     alpha = alpha, 
                                     pfs = pf[feature.cut], 
                                     lambda.min.ratio = lambda.min.ratio, 
                                     nlam = nlam,
                                     weights = weights)
    }else{
      lambda = get.lambda(X[,feature.cut], y, 
                          family = family,
                          alpha = alpha,
                          nlam = nlam, 
                          lambda.min.ratio = lambda.min.ratio,
                          weights = weights,
                          pfs = pf[feature.cut])
    }
  }
  
  fit <- glmnet(x = X[,feature.cut], y = y,
                family = family,
                alpha = alpha,
                lambda = lambda, 
                penalty.factor = pf[feature.cut],
                weights = weights,
                ...)
  fit$cut = cut
  fit$feature.cut = feature.cut
  fit
}

#-------------------------------------------------------------------------------


