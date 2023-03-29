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
#' We take in:
#'  1. List of formulae (list of lists, with elements fixed, random, response)
#'     These are NOT in formula form.
#'  2. The data
#' So, let's try::
#' 
library(Rcpp); library(RcppArmadillo)
sourceCpp('dev/cpp-modelmatrices.cpp')

AssignIds <- function(data){
  id <- data$id
  uid <- unique(data$id)
  list(id = uid, assign = 1:length(uid))
}

SubsetData <- SubsetAllData(data, AssignIds(data))

long.formulas <- list(
  Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
  Y.2 ~ time + cont + bin + (1 + time|id)  # Poisson
)

formulas <- lapply(long.formulas, parseFormula)

CreateDataMatricesNew <- function(formulas, data){
  # Data subsets
  idList <- AssignIds(data)
  idata <- SubsetAllData(data, idList)
  n <- length(idList$id)
  K <- length(formulas)
  
  # Transform formulae made by parseFormula into `~ x` formula objects
  formulas2 <- lapply(formulas, function(f){
    f$fixed <- structure(as.formula(paste0('~ ', f$fixed)), special = attr(f$fixed, 'special'))
    f$random <- structure(as.formula(paste0('~ ', f$random)), special = attr(f$random, 'special'))
    f
  })
  
  D <- lapply(formulas2, function(f){
    NoSplineCheck <- attr(f$fixed, 'special') == 'none'
    if(NoSplineCheck){
      X <- Map(function(idat) CreateMatrices(f$fixed, idat, "X", f$resp), idat = idata)
      Z <- Map(function(idat) CreateMatrices(f$random, idat, "Z", f$resp), idat = idata)
      Y <- Map(function(idat) CreateMatrices(f$fixed, idat, "Y", f$resp), idat = idata)
    }else{
      newData <- as.data.frame(cbind(id = data$id, model.matrix(f$fixed, data)))
      newDataRandom <- as.data.frame(cbind(id = data$id , model.matrix(f$random, data)))
      X <- Z <- vector("list", length = length(idList$assign))
      for(i in seq_along(idList$id)){
        X[[idList$assign[i]]] <- newData[newData$id == idList$id[i], ]
        Z[[idList$assign[i]]] <- newDataRandom[newDataRandom$id == idList$id[i], ]
      }
    }
    Y <- Map(function(idat) CreateMatrices(f$fixed, idat, "Y", f$resp), idat = idata)
    list("X" = X, "Z" = Z, "Y" = Y)
  })
  # Creates objects of class "AsIs" if e.g. quadratic term is supplied, 
  # Not sure on ramifications here.
  
  X <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      as.matrix(D[[k]]$X[[i]])
    })
  })
  Z <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      as.matrix(D[[k]]$Z[[i]])
    })
  })
  Y <- lapply(1:n, function(i){
    lapply(1:K, function(k){
      as.matrix(D[[k]]$Y[[i]])
    })
  })
  
  P <- sapply(X[[1]], ncol)
  nP <- lapply(X[[1]], colnames)
  q <- sapply(Z[[1]], ncol)
  nq <- lapply(Z[[1]], colnames)
  
  list(X = X, Y = Y, Z = Z,
       P = P, nP = nP, q = q, nq = nq)
}

old.vs.new <- microbenchmark::microbenchmark(
  `old` = {
    createDataMatrices(data, formulas)
  },
  `new` = {
    CreateDataMatricesNew(formulas, data)
  },
  times = 50L
)

# Check same
old <- createDataMatrices(data, formulas)
new <- CreateDataMatricesNew(formulas, data)

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

old <- createDataMatrices(data, formulas10)
new <- CreateDataMatricesNew(formulas10, data)

old.vs.new10 <- microbenchmark::microbenchmark(
  `old` = {
    createDataMatrices(data, formulas10)
  },
  `new` = {
    CreateDataMatricesNew(formulas10, data)
  },
  times = 50L
)
t.test(old.vs.new10$time[old.vs.new10$expr=='old'], old.vs.new10$time[old.vs.new10$expr=='new'])