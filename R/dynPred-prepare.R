# Preparing data for dynamic prediction routine.

#' @keywords internal
prepareLongData <- function(data, long.formulas, responses, u = NULL){
  # Prepares data for one subject i. That is, data argument should be data[data$id==i, ] type
  # long.formulas/responses should come from fitted `joint` object. (i.e. joint$modelInfo)
  # If u is NULL then the whole data[,] object is partitioned and returned.
  #     One should use the `u` argument to project from Tstart in the Tstart+delta window.
  
  formulas <- lapply(long.formulas, parseFormula); 
  K <- length(formulas)
  responses <- lapply(strsplit(responses, '\\s\\('), el, 1)
  
  if(!is.null(u)) data <- data[data$time <= u, ] # Truncate if necessary
  
  X <- lapply(1:K, function(k){
    .DataWorkhorse(data, data$id[1], formulas[[k]]$fixed, formulas[[k]]$random, responses[[k]], what = 'X')
  })
  Y <- lapply(1:K, function(k){
    .DataWorkhorse(data, data$id[1], formulas[[k]]$fixed, formulas[[k]]$random, responses[[k]], what = 'Y')
  })
  Z <- lapply(1:K, function(k){
    .DataWorkhorse(data, data$id[1], formulas[[k]]$fixed, formulas[[k]]$random, responses[[k]], what = 'Z')
  })
  
  list(X = X, Y = Y, Z = Z)
}

#' @keywords internal
prepareSurvData <- function(data, fit, formulas, u = NULL, hazard){
  # data: i.dat-type object; survival data is then derived in-formula (by taking the first row).
  # fit: joint model fit
  # formulas: list of formulae from the fitted aEM object.
  # u: horizon time; if not supplied then this is found for the entire profile is returned.
  # hazard: hazard from fitted aEM object.
  
  formulas <- lapply(formulas, parseFormula); K <- length(formulas)
  
  # Extract failure times and baseline hazard
  ft <- hazard[,1]; l0 <- hazard[,2]
  
  survtime <- data$survtime[1]
  status <- as.integer(survtime <= max(data$time)) # Do they die in the window 0 < t / 0 < u?
  data$longtime <- data$time; max.longtime <- max(data$longtime)
  if(max.longtime == 0) max.longtime <- data$survtime # prevent issues
  data$time <- data$survtime # potentially bad move, we'll see!
  data <- data[1, ]
 
  if(!is.null(u)) data$time <- u # Don't extrapolate beyond time u
  if(is.null(u) && data$time > max.longtime) data$time <- max.longtime
  
  # Fi
  Fi <- do.call(cbind, lapply(formulas, function(x){
    model.matrix(as.formula(paste0('~', x$random)), data)
  }))
  
  # Fu and l0u
  # Work out survived times (and if necessary truncate at u if supplied)
  survived.times <- ft[which(ft <= data$time)]
  
  if(length(survived.times) > 0){ # Avoid complications if they're censored before first failure time.
    Fu <- do.call(cbind, lapply(formulas, function(x){
      model.matrix(as.formula(paste0('~', x$random)), data.frame(time=survived.times))
    }))
    l0u <- l0[1:length(survived.times)]
  }else{
    Fu <- matrix(0, nr = 1, nc = ncol(Fi))
    l0u <- 0
  }
  
  if(status == 1L){
    l0i <- l0[which(ft == survtime)]
    if(!is.null(u) && u < survtime) l0i <- l0[max(which(ft <= u))] # truncate if necessary
  }else{
    l0i <- 0
  }
  Sname <- gsub("^zeta\\_", '', names(fit$coeffs$zeta))
  S <- as.matrix(data[1,Sname,drop=F])
  SS <- apply(S, 2, rep, nrow(Fu))
  
  list(
    S = S, SS = SS,
    l0u = l0u, l0i = l0i,
    Fi = Fi, Fu = Fu, Delta = status
  )
}

# Prepare and collate longitudinal and survival data for dynamic predictions
#' @keywords internal
prepareData <- function(data, id, fit, u = NULL){
  long.formulas <- fit$ModelInfo$long.formulas
  response.info <- fit$ModelInfo$ResponseInfo
  long <- prepareLongData(data, long.formulas, response.info, u = u)
  surv <- prepareSurvData(data, fit, long.formulas, u = u, hazard = fit$hazard)
  
  b <- fit$REs[id,]; K <- length(fit$family)
  responsenames <- lapply(strsplit(response.info, '\\s\\('), el , 1)
  
  b.inds <- fit$ModelInfo$inds$b
  gamma.rep <- rep(fit$coeffs$gamma, sapply(b.inds, length))
  
  # Ensure zero-indexing in call to joint_density.
  b.inds <- lapply(b.inds, function(x) x - 1)     
  beta.inds <- lapply(fit$ModelInfo$inds$beta, function(x) x - 1)
  
  if(is.null(u)){
    bfit <- optim(
      rep(0, ncol(fit$REs)), joint_density, joint_density_ddb, Y = long$Y, X = long$X, Z = long$Z,
      beta = fit$coeffs$beta, D = fit$coeffs$D, sigma = fit$coeffs$sigma, family = fit$family,
      Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i, SS = surv$SS, Fu = surv$Fu,
      haz = surv$l0u, gamma_rep = gamma.rep, zeta = fit$coeffs$zeta, beta_inds = beta.inds,
      b_inds = b.inds, K = K, method = 'BFGS', hessian = T
    )
  }else{
    bfit <- NULL
  }
  
  list(
    long = long, surv = surv, b.MLE = b, b = bfit
  )
  
}