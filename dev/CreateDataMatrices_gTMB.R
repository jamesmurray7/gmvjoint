# R/CreateDataMatrices.R --------------------------------------------------
.DataWorkhorse <- function(data, subj.id, fixed, random, response, what = 'X'){
  fixed.formula <- as.formula(paste0('~', fixed))
  random.formula <- as.formula(paste0('~', random))
  if(what == 'Y') rtn <- matrix(data[data$id == subj.id, response], ncol = 1)
  else{
    if(!is.null(attr(fixed, 'special')) & (attr(fixed, 'special') == 'spline' | attr(random, 'special') == 'spline')){
      newData <- as.data.frame(cbind(id = data$id, model.matrix(fixed.formula, data)))
      newDataRandom <- as.data.frame(cbind(id = data$id, model.matrix(random.formula, data)))
      if(what == 'X') rtn <- as.matrix(newData[newData$id == subj.id, -1, drop = F])
      else if(what == 'Z') rtn <- as.matrix(newDataRandom[newDataRandom$id == subj.id, -1, drop = F])
    }else{ 
      i.dat <- data[data$id == subj.id, ]
      if(what == 'X') rtn <- model.matrix(fixed.formula, i.dat)
      else if(what == 'Z') rtn <- model.matrix(random.formula, i.dat)
    }
  }
  rtn
}

createDataMatrices <- function(data, formulas){
  K <- length(formulas)
  n <- length(unique(data$id))
  X <- Y <- Z <- vector('list', n)
  for(i in 1:n){
    ii <- i
    X[[i]] <- lapply(formulas, function(f){
      .DataWorkhorse(data, subj.id = ii, fixed = f$fixed, random = f$random, response = f$response, what = 'X')
    }) 
    Z[[i]] <- lapply(formulas, function(f){
      .DataWorkhorse(data, subj.id = ii, fixed = f$fixed, random = f$random, response = f$response, what = 'Z')
    }) 
    Y[[i]] <- lapply(formulas, function(f){
      .DataWorkhorse(data, subj.id = ii, fixed = f$fixed, random = f$random, response = f$response, what = 'Y')
    }) 
  }
  
  list(X = X, Y = Y, Z = Z,
       P = sapply(X[[1]], ncol),
       q = sapply(Z[[1]], ncol),
       np = lapply(X[[1]], colnames),
       nq = lapply(Z[[1]], colnames)
  )
}


# Setting out goals -------------------------------------------------------
#' We obtain a model fit by glmmTMB as part of model wind-up. 
#' How easily can we form requisite data matrices?
#' And is this faster than the old method (given above)?

AssignIds <- function(data){
  id <- data$id
  uid <- unique(data$id)
  list(id = uid, 
       num.each = unname(table(id)),
       assign = 1:length(uid),
       n = length(uid), rawid = id)
}


long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
  Y.2 ~ time + cont + bin + (1 + time|id)   # Poisson
)

attr(long.formulas[[1]], "family") <- gaussian()
attr(long.formulas[[2]], "family") <- poisson()

fits <- lapply(long.formulas, function(f){
  glmmTMB(f, data = data,
          family = attr(f, 'family'))
})

idList <- AssignIds(data)
.getdmats <- function(f, idList){
  stopifnot("glmmTMB" %in% class(f))
  q.names <- f$modelInfo$reTrms$cond$cnms$id
  resp <- names(f$modelInfo$respCol)
  q <- length(q.names)
  allX <- getME(f, "X")
  allZ <- as.matrix(getME(f, "Z"))
  allY <- matrix(f$frame[, resp], ncol = 1)
  colnames(allY) <- resp
  inds <- split(1:nrow(allX), idList$rawid)
  col.inds <- split(colnames(allZ), rep(idList$id, each = q))
  col.inds.w <- lapply(col.inds, function(x) which(colnames(allZ)==x))
  colnames(allZ) <- rep(q.names, idList$n)
  X <- lapply(inds, function(i) allX[i,])
  Z <- lapply(seq_along(inds), function(i) allZ[inds[[i]], col.inds.w[[i]]])
  Y <- lapply(inds, function(i) allY[i,])
  list(X = X, 
       Y = Y, 
       Z = Z,
       q = q,
       q.names = q.names,
       P = ncol(allX),
       P.names = colnames(allX))
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
  
  P <- sapply(1:K, function(k) D[[k]]$P)
  q <- sapply(1:K, function(k) D[[k]]$q)
  nP <- lapply(1:K, function(k) D[[k]]$P.names)
  nq <- lapply(1:K, function(k) D[[k]]$q.names)
  
  list(X = X, Y = Y, Z = Z,
       P = P, nP = nP, q = q, nq = nq)
}

dmats <- getdmats(fits, idList)

formulas <- lapply(long.formulas, parseFormula)
old.vs.new <- microbenchmark::microbenchmark(
  `old` = {
    createDataMatrices(data, formulas)
  },
  `new` = {
    dmats <- getdmats(fits, idList)
  },
  times = 50L
)

# Check same
old <- createDataMatrices(data, formulas)
new <- getdmats(fits, idList)

old.vs.new
plot(old.vs.new)

# What about a 10-variate?
data <- simData(beta = do.call(rbind, replicate(10, c(2,-1,0.1, .1), simplify = F)),
                D = diag(rep(c(.25, .04), 10)),
                gamma = sample(c(.5, -.5), 10, T),
                family = as.list(rep("gaussian", 10)),
                sigma = rep(.16, 10))$data

long.formulas10 <- lapply(1:10, function(k){
  as.formula(paste0('Y.', k, ' ~ time + cont + bin + (1 + time|id)'))
})
formulas10 <- lapply(long.formulas10, parseFormula)

for(k in 1:10) attr(long.formulas10[[k]], 'family') <- 'gaussian'

fits10 <- lapply(long.formulas10, function(f){
  glmmTMB(f, data = data,
          family = attr(f, 'family'))
})

idList10 <- AssignIds(data)
old <- createDataMatrices(data, formulas10)
new <- getdmats(fits10, idList10)

old.vs.new10 <- microbenchmark::microbenchmark(
  `old` = {
    old <- createDataMatrices(data, formulas10)
  },
  `new` = {
    dmats <- getdmats(fits10, idList10)
  },
  times = 50L
)
old.vs.new10 # Very large timesave!
plot(old.vs.new10)
t.test(old.vs.new10$time[old.vs.new10$expr=='old'], old.vs.new10$time[old.vs.new10$expr=='new'])
