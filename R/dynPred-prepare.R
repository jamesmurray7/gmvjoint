# Preparing data for dynamic prediction routine.

#' @keywords internal
prepareLongData <- function(data, fit){
  
  K <- fit$ModelInfo$K
  id <- unique(data$id)
  
  X <- fit$dmats$long$X[[id]]
  Z <- fit$dmats$long$Z[[id]]
  Y <- fit$dmats$long$Y[[id]]
  W <- fit$dmats$long$W[[id]]
  
  # Truncate at max. time value in `data`, which should look like
  # data[data$time == x, ].
  
  # trunc.time <- max(data$time)
  ln <- nrow(data)
  
  X <- lapply(1:K, function(k){
    Xk <- X[[k]]
    # Xk[Xk[, 'time'] <= trunc.time, , drop = F]
    Xk[1:ln, , drop=F]
  })
  Z <- lapply(1:K, function(k){
    Zk <- Z[[k]]
    # Zk[Zk[, 'time'] <= trunc.time, , drop = F]
    Zk[1:ln, , drop = F]
  })
  Y <- lapply(1:K, function(k){
    # Y[[k]][1:nrow(X[[k]])]
    Y[[k]][1:ln]
  })
  W <- lapply(1:K, function(k){
    Wk <- W[[k]]
    # Wk[1:nrow(X[[k]]),,drop=F]
    Wk[1:ln,,drop=F]
  })
  
  list(X = X, Y = Y, Z = Z, W = W)
}

#' @keywords internal
prepareSurvData <- function(data, fit, u = NULL){
  
  # This id 
  id <- unique(data$id)
  if(is.null(u)){ # Set to maximum observed longitudinal time if not defined.
    u <- max(data$time)
  }
  
  # All _possible_ survived times observed in model, cast by W_k(t).
  # along with baseline hazard, as appeared in original model...
  Fu.all <- fit$dmats$surv$ft.mat
  l0u <- fit$hazard[,2]
  ft <- fit$dmats$surv$ft
  
  fs <- lapply(fit$ModelInfo$long.formulas, parseFormula)
  Wk <- lapply(fs, '[[', 'random')
  
  # Need to re-derive for each u
  Fi <- do.call(cbind, lapply(seq_along(Wk), function(k){
    spec <- attr(Wk[[k]], 'special')
    if(spec == "none"){
      out <- model.matrix(as.formula(paste0("~", Wk[[k]])), data.frame(time = u))
    }else{
      sft <- fit$dmats$surv$spline.fts[[k]]
      out <- cbind(1, predict(sft, data.frame(time = u)))
    }
    out
  }))
  
  Fu <- Fu.all[ft <= u,]
  if(nrow(Fu)){ # Code-in a failsafe.
    haz <- l0u[1:nrow(Fu)]
  }else{
    Fu <- matrix(0, nrow = 1, ncol = ncol(Fu.all))
    haz <- c(0)
  }
  
  # Fail safe for Fu
  
  # Time-invariant items
  S <- fit$dmats$surv$S[[id]]
  SS <- apply(S, 2, rep, nrow(Fu))
  if(!"matrix"%in%class(SS)) SS <- as.matrix(t(SS))
  
  return(list(
    S = S, SS = SS, 
    Fi = Fi, Fu = Fu, l0u = haz,
    l0i = 0., Delta = 0L     # Suppose we don't know T_i*
  ))
  
}


# Prepare and collate longitudinal and survival data for dynamic predictions
#' @keywords internal
prepareData <- function(data, fit, u = NULL){
  id <- unique(data$id)
  
  if(is.null(u)) long <- prepareLongData(data, fit) # Only needed in denominator set.
  surv <- prepareSurvData(data, fit, u = u)
  
  if(is.null(u)){ # This to work out b.current and Sigma at start of MH.
    b <- fit$REs[id, ]
    K <-fit$ModelInfo$K
    b.inds <- fit$ModelInfo$inds$C$b; gamma.rep <- rep(fit$coeffs$gamma, sapply(b.inds, length))
    beta.inds <- fit$ModelInfo$inds$C$beta
    
    bfit <- optim(rep(0, ncol(fit$REs)),
                  joint_density,
                  joint_density_ddb,
                  Y = long$Y, X = long$X, Z = long$Z, W = long$W, beta = fit$coeffs$beta,
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
