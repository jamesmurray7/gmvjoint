# Quite slow in higher dimension; could work out a C++ alternative?
#' @keywords internal
.DataWorkhorse <- function(data, subj.id, fixed, random, response, what = 'X'){
  fixed.formula <- as.formula(paste0('~', fixed))
  random.formula <- as.formula(paste0('~', random))
  if(what == 'Y') rtn <- matrix(data[data$id == subj.id, response],nc = 1)
  else{
    if(!is.null(attr(fixed, 'special')) & (attr(fixed, 'special') == 'spline' | attr(random, 'special') == 'spline')){
      newData <- as.data.frame(cbind(id = data$id, model.matrix(fixed.formula, data)))
      newDataRandom <- as.data.frame(cbind(id = data$id, model.matrix(random.formula, data)))
      if(what == 'X') rtn <- as.matrix(newData[newData$id == subj.id, -1, drop = F])
      else if(what == 'Z') rtn <- as.matrix(newDataRandom[newDataRandom$id == subj.id, -1, drop = F])
    }else{ 
      i.dat <- subset(data, id == subj.id)
      if(what == 'X') rtn <- model.matrix(fixed.formula, i.dat)
      else if(what == 'Z') rtn <- model.matrix(random.formula, i.dat)
    }
  }
  rtn
}

# Create lists of design matrices for each subject for each response
#' @keywords internal
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