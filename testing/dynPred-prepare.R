# TESTING VERSION -- with dmats
# Preparing data for dynamic prediction routine.

# Prepare longitudinal data for dynamic predictions.
# This will only _ever_ be used in denominator of S(u|t)/S(t|t)
# where t is defined by first element of u vector in dynPreds.
prepareLongData <- function(data, fit){
  
  K <- length(fit$ModelInfo$long.formulas)
  id <- unique(data$id)
  
  X <- fit$dmats$long$X[[id]]
  Z <- fit$dmats$long$Z[[id]]
  Y <- fit$dmats$long$Y[[id]]
  
  # Truncate at max. time value in `data`, which should look like
  # data[data$time == x, ].
  
  trunc.time <- max(data$time)
  X <- lapply(1:K, function(k){
    Xk <- X[[k]]
    Xk[Xk[, 'time'] <= trunc.time, ]
  })
  Z <- lapply(1:K, function(k){
    Zk <- Z[[k]]
    Zk[Zk[, 'time'] <= trunc.time, ]
  })
  Y <- lapply(1:K, function(k){
    Y[[k]][1:nrow(X[[k]])]
  })
  
  list(X = X, Y = Y, Z = Z)
}

prepareSurvData <- function(data, fit, u = NULL){
  
  # This id 
  id <- unique(data$id)
  if(is.null(u)){ # Set to maximum observed longitudinal time if not defined.
    u <- max(data$time)
  }
  
  # All _possible_ survived times observed in model, cast by W_k(t).
  # along with baseline hazard
  Fu.all <- fit$dmats$surv$Fu[[which.max(sapply(fit$dmats$surv$Fu, nrow))]]
  l0u <- fit$hazard[,2]
  ft <- fit$dmats$surv$ft
  
  fs <- lapply(fit$ModelInfo$long.formulas, parseFormula)
  Fi <- do.call(cbind, lapply(fs, function(x) model.matrix(as.formula(paste0('~', unlist(x$random))),
                                      data.frame(time = u)))) # Don't actually need to work this out.
  Fu <- Fu.all[ft <= u,]; haz <- l0u[1:nrow(Fu)]
  
  # Time-invariant items
  S <- fit$dmats$surv$S[[id]]
  SS <- fit$dmats$surv$SS[[id]][1:nrow(Fu), , drop = FALSE]
  
  return(list(
    S = S, SS = SS, 
    Fi = Fi, Fu = Fu, l0u = haz,
    l0i = 0., Delta = 0L     # Supppose we don't know T_i*
  ))
  
}

prepareData <- function(data, fit, u = NULL){
  id <- unique(data$id)
  
  if(is.null(u)) long <- prepareLongData(data, fit) # Only needed in denominator set.
  surv <- prepareSurvData(data, fit, u = u)
  
  if(is.null(u)){ # This to work out b.current and Sigma at start of MH.
    b <- fit$REs[id, ]
    K <- length(fit$ModelInfo$family)
    b.inds <- fit$ModelInfo$inds$b; gamma.rep <- rep(fit$coeffs$gamma, sapply(b.inds, length))
    b.inds <- lapply(b.inds, function(x) x - 1)
    beta.inds <- lapply(fit$ModelInfo$inds$beta, function(x) x - 1)
    
    bfit <- optim(rep(0, ncol(fit$REs)),
                  joint_density,
                  joint_density_ddb,
                  Y = long$Y, X = long$X, Z = long$Z, beta = fit$coeffs$beta,
                  D = fit$coeffs$D, sigma = fit$coeffs$sigma, family = fit$ModelInfo$family,
                  Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
                  SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma.rep,
                  zeta = fit$coeffs$zeta, beta_inds = beta.inds, b_inds = b.inds, K = K,
                  method = 'BFGS', hessian = TRUE)
    
    if(bfit$convergence!=0L) 
      warning(sprintf("Unsuccessful minimisation of longitudinal data trucated at chosen value (convergence code: %d)",
                      bfit$convergence))
  }
  
  out <- list()
  out$surv <- surv
  if(is.null(u)){
    out$long <- long
    out$b.MLE <- b
    out$bfit <- bfit
  }
  
  out
}
