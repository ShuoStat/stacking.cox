
#' cv.drpcr 
#' 
#' @description 
#' 
#' Cross-validation function for dimension reduction based on 
#' principle components, a method to combine clinical 
#' molecular variables for prediction. Step1, derive principle components; 
#' Step2, model fitting.
#' 
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' @param blocks block of variables
#' @param blocks.pc block variables for principle components
#' @param foldid list or vector depends on ncv
#' @param nfolds defaulted as 10
#' @param alpha defaulted as 1 
#' @param penalty.factors list of penalty factors 
#' @param type.measure "default" 
#' @param lambda NULL 
#' @param nlam 100 
#' @param lambda.min.ratio 0.001 
#' 
#' @examples 
#' 
#' @import glmnet
#' @export cv.drpcr
#' 

cv.drpcr <- function(X, y, alpha = 1, family,
                     ncomp = nrow(X) - 1,
                     blocks = rep(1, ncol(X)), 
                     blocks.pc = c(1, 2),
                     weights = rep(1, nrow(X)),
                     penalty.factors = NULL,
                     foldid = NULL, 
                     nfolds = 10, 
                     lambda = NULL, 
                     lambda.min.ratio = 0.001, 
                     nlam = 100, 
                     type.measure = "default"){
  
  n = nrow(X)
  p = ncol(X)
  
  nblocks = length(unique(blocks))
  #- get blocks
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
  
  #- blocks
  #- centers or sds
  cts <- rep(0, length(blocks)) #- centers
  sds <- rep(1, length(blocks)) #- sds
  new.blocks <- c()
  M <- c() #- new data
  for(i in unique(blocks)){
    
    s <- blocks == i
    Z <- X[, s]
    if (i %in% blocks.pc) {
      
      #- pcr, normalization and scaling
      Z <- scale(Z, scale = T, center = T)
      cts[s] <- attr(Z, "scaled:center")
      sds[s] <- attr(Z, "scaled:scale")
      SVD <- svd(Z)
      tmp <- Z %*% SVD$v[,1:ncomp]
      colnames(tmp) <- paste0("pc", i, "_", seq_len(ncol(tmp)))
      
    } else {
      tmp <- Z
    }
    M <- cbind(M, tmp)
    new.blocks <- c(new.blocks, rep(i, ncol(tmp)))
  }
  
  if (is.null(foldid))
    foldid <- sample(rep(1:nfolds, length = n), n)
  
  #- using cv.ipflas in fav package 
  cvs <- cv.ipflas(M, y, 
                   family = family, 
                   foldid = foldid, 
                   nfolds = nfolds, 
                   alpha = alpha, 
                   weights = weights,
                   blocks = new.blocks,
                   penalty.factors = penalty.factors, 
                   type.measure = type.measure,
                   lambda = lambda, 
                   nlam = nlam, 
                   lambda.min.ratio = lambda.min.ratio)
  
  append(cvs, list(centers = cts, scales = sds, blocks = new.blocks,
                   ncomp = ncomp))
  
}


#' drpcr 
#' 
#' @description 
#' 
#' Function for dimension reduction based on principle components
#' 
#' @param X matrix
#' @param y vector of response
#' @param family defaulted as gaussian
#' @param blocks block of variables
#' @param blocks.pc block variables for principle components
#' @param nfolds defaulted as 10
#' @param alpha defaulted as 1 
#' @param pf penalty factor
#' @param type.measure "default" 
#' @param lambda NULL 
#' @param nlam 100 
#' @param lambda.min.ratio 0.001 
#' 
#' @examples 
#' 
#' @import glmnet
#' @import Matrix
#' @export drpcr
#' 
#' 
#' 
drpcr <- function(X, y, family,
                  ncomp = NULL,
                  alpha = 1, 
                  blocks = rep(1, ncol(X)), 
                  blocks.pc = c(1),
                  weights = rep(1, nrow(X)),
                  pf = NULL,
                  lambda = NULL, lambda.min.ratio = 0.001, nlam = 100, 
                  type.measure = "default"){
  
  #- blocks
  #- centers or sds
  cts <- rep(0, length(blocks)) #- centers
  sds <- rep(1, length(blocks)) #- sds
  new.blocks <- c()
  
  M <- c() #- new data
  eigin.vectors <- list()
  for(i in unique(blocks)){
    
    s <- blocks == i
    Z <- X[, s]
    if (i %in% blocks.pc) {
      #- pcr, normalization and scaling
      Z <- scale(Z, scale = T, center = T)
      cts[s] <- attr(Z, "scaled:center")
      sds[s] <- attr(Z, "scaled:scale")
      SVD <- svd(Z)
      tmp <- Z %*% SVD$v[,1:ncomp]
      colnames(tmp) <- paste0("pc", i, "_", seq_len(ncol(tmp)))
      eigin.vectors[[i]] <- SVD$v[,1:ncomp]
    } else {
      tmp <- Z
    }
    M <- cbind(M, tmp)
    new.blocks <- c(new.blocks, rep(i, ncol(tmp)))
  }
  
  ipf <- ipflas(M, y, 
                family = family, 
                alpha = alpha, 
                blocks = new.blocks, 
                pf = pf, 
                lambda = lambda,
                weights = weights,
                nlam = nlam, 
                lambda.min.ratio = lambda.min.ratio)
  
  beta <- ipf$beta
  beta0 <- Matrix::Matrix(nrow = ncol(X), ncol = ncol(beta),
                          data = 0, sparse = T)
  
  for(i in unique(blocks)){
    
    s <- new.blocks == i
    if(i %in% blocks.pc){
      v <- eigin.vectors[[i]]
      tmp <- v %*% beta[s,] 
    }else{
      tmp <- beta[s,]
    }
    beta0[blocks == i,] <- tmp
  }
  
  beta0 <- beta0 / sds
  a0 <- ipf$a0 - t(beta0) %*% cts
  
  back.mod <- list(a0 = a0, beta = beta0)
  out <- append(ipf, list(centers = cts, 
                          scales = sds, 
                          blocks = new.blocks,
                          variables = colnames(X),
                          back.mod = back.mod))
  class(out) <- "drpcr"
  out
}

#-------------------------------------------------------------------------------

drpcr2 <- function(X, y, alpha = 1, family,
                   ncomp = ncomp,
                   blocks = rep(1, ncol(X)), 
                   blocks.pc = c(1),
                   weights = weights,
                   pf = NULL,
                   lambda = NULL, lambda.min.ratio = 0.001, nlam = 100, 
                   type.measure = "default"){
  
  #- blocks
  #- centers or sds
  cts <- rep(0, length(blocks)) #- centers
  sds <- rep(1, length(blocks)) #- sds
  new.blocks <- c()
  
  M <- c() #- new data
  eigin.vectors <- list()
  for(i in unique(blocks)){
    
    s <- blocks == i
    Z <- X[, s]
    if (i %in% blocks.pc) {
      #- pcr, normalization and scaling
      Z <- scale(Z, scale = T, center = T)
      cts[s] <- attr(Z, "scaled:center")
      sds[s] <- attr(Z, "scaled:scale")
      SVD <- svd(Z)
      tmp <- Z %*% SVD$v[,1:ncomp]
      colnames(tmp) <- paste0("pc", i, "_", seq_len(ncol(tmp)))
      eigin.vectors[[i]] <- SVD$v[,1:ncomp]
    } else {
      tmp <- Z
    }
    M <- cbind(M, tmp)
    new.blocks <- c(new.blocks, rep(i, ncol(tmp)))
  }
  
  ipf <- ipflas(M, y, 
                family = family, 
                alpha = alpha, 
                blocks = new.blocks, 
                pf = pf, 
                lambda = lambda,
                weights = weights, 
                nlam = nlam, 
                lambda.min.ratio = lambda.min.ratio)
  
  out <- append(ipf, list(centers = cts, 
                          scales = sds, 
                          blocks = new.blocks,
                          variables = colnames(X),
                          blocks.pc = blocks.pc,
                          eigin.vectors = eigin.vectors))
  class(out) <- "drpcr"
  out
}

predict.drpcr <- function(obj, newX, 
                          type = c("link", "response", "coefficients", 
                                   "nonzero", "class")){
  
  eigin.vectors = obj$eigin.vectors
  centers <- obj$centers
  scales <- obj$scales
  M <- c() #- new data
  for(i in unique(blocks)){
    
    s <- blocks == i
    Z <- X[, s]
    if (i %in% blocks.pc) {
      
      #- pcr, normalization and scaling
      Z <- scale(Z, center = centers[s], scale = scales[s])
      ens <- eigin.vectors[[i]]
      tmp <- Z %*% ens
      colnames(tmp) <- paste0("pc", i, "_", seq_len(ncol(tmp)))
    } else {
      tmp <- Z
    }
    M <- cbind(M, tmp)
  }
  
  predict.glmnet(obj, M, type = type)
}

#- end -------------------------------------------------------------------------















