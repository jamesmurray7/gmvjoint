AssignIds <- function(data){
  id <- data$id
  uid <- unique(data$id)
  list(id = uid, 
       num.each = unname(table(id)),
       assign = 1:length(uid),
       n = length(uid), rawid = id)
}

.getdmats <- function(f, idList){
  stopifnot("glmmTMB" %in% class(f))
  q.names <- f$modelInfo$reTrms$cond$cnms$id
  resp <- names(f$modelInfo$respCol)
  q <- length(q.names)
  allX <- glmmTMB::getME(f, "X")
  allW <- glmmTMB::getME(f, "Xd")
  allZ <- as.matrix(glmmTMB::getME(f, "Z"))
  allY <- matrix(f$frame[, resp], ncol = 1)
  colnames(allY) <- resp
  inds <- split(1:nrow(allX), idList$rawid)
  col.inds <- split(colnames(allZ), rep(idList$id, each = q))
  col.inds.w <- lapply(col.inds, function(x) which(colnames(allZ)==x))
  colnames(allZ) <- rep(q.names, idList$n)
  # {X, W, Z} should be matrices, so enforce drop = F.
  X <- lapply(inds, function(i) allX[i,,drop=F])
  W <- lapply(inds, function(i) allW[i,,drop=F])
  Z <- lapply(seq_along(inds), function(i) allZ[inds[[i]], col.inds.w[[i]], drop = F])
  Y <- lapply(inds, function(i) allY[i,])
  list(X = X, 
       W = W,
       Y = Y, 
       Z = Z,
       q = q,
       q.names = q.names,
       P = ncol(allX),
       Pd = ncol(allW),
       P.names = colnames(allX),
       Pd.names = colnames(allW))
}

getdmats <- function(fits, idList){
  D <- lapply(fits, .getdmats, idList)
  
  K <- length(fits)
  n <- idList$n
  # Rearrange
  X <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      D[[k]]$X[[i]]
    })
  })
  W <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      D[[k]]$W[[i]]
    })
  })
  Y <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      D[[k]]$Y[[i]]
    })
  })
  Z <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      D[[k]]$Z[[i]]
    })
  })
  # Total num. obs per response.
  m <- sapply(1:K, function(k) fits[[k]]$modelInfo$nobs)
  # Number of observations per subject per response.
  mi <- lapply(Y, function(y) lapply(1:K, function(k) length(y[[k]])))
  
  P <- sapply(1:K, function(k) D[[k]]$P)
  Pd <- sapply(1:K, function(k) D[[k]]$Pd)
  q <- sapply(1:K, function(k) D[[k]]$q)
  nP <- lapply(1:K, function(k) D[[k]]$P.names)
  nPd <- lapply(1:K, function(k) D[[k]]$Pd.names)
  nq <- lapply(1:K, function(k) D[[k]]$q.names)
  
  list(X = X, Y = Y, Z = Z, W = W,
       P = P, nP = nP, 
       Pd = Pd, nPd = nPd,
       q = q, nq = nq,
       K = K, n = n, m = m, mi = mi)
}