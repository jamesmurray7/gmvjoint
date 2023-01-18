#' @keywords internal
#' @importFrom MASS mvrnorm
#' @importFrom Matrix nearPD
#' @importFrom stats vcov
Omega.draw <- function(x){
  Omega.Var <- vcov(x)
  # Mean vector
  co <- x$coeffs
  sigma <- unlist(co$sig)
  sigma <- sigma[sigma != 0.0]
  Omega.mean <- setNames(
    c(
      c(vech(co$D)),
      co$beta, sigma,
      co$gamma, co$zeta
    ), names(x$SE))
  
  # Draw from N(Omega.mean, Omega.Var)
  draw <- MASS::mvrnorm(n = 1, mu = Omega.mean, Sigma = Omega.Var)
  
  # Re-construct Omega at this draw.
  D <- matrix(0, nrow(co$D), ncol(co$D)) # Need to check this works for K>1.
  D[lower.tri(D, T)] <- draw[grepl('^D\\[', names(draw))]
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  # Check this is pos-definite and transform if not.
  if(any(eigen(D)$values < 0) || (det(D) <= 0)){ 
    D <- as.matrix(Matrix::nearPD(D)$mat)
  }
  beta <- draw[match(names(co$beta), names(draw))]
  gamma <- draw[match(names(co$gamma), names(draw))]
  zeta <- draw[match(names(co$zeta), names(draw))]
  # Sort out sigma --> MUST be a list.
  .sigma <- lapply(co$sigma, function(x){
    if(x == 0L)
      return(0)
    else
      return(draw[match(names(x), names(draw))])
  })
  
  list(
    D = D, beta = beta, sigma = .sigma, gamma = gamma, zeta = zeta
  )
}

# Complete data log-likelihood w.r.t random effects b.
#' @keywords internal
logLik.b <- function(b, long, surv, O, beta.inds, b.inds, fit){
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(b.inds, length)); zeta <- O$zeta; sigma <- O$sigma
  
  neg.ll <- joint_density(b = b, Y = long$Y, X = long$X, Z = long$Z, beta = beta, D = D, sigma = sigma,
                          family = fit$ModelInfo$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
                          SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
                          b_inds = b.inds, K = length(fit$ModelInfo$family))
  
  -neg.ll
}

# Metropolis-Hastings scheme to draw random effects from candidate density.
#' @keywords internal
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
b.mh <- function(b.current, b.hat.t, Sigma.t, long, surv, O, beta.inds, b.inds, fit, 
                 b.density, df){
  # Metropolis-Hasting scheme to draw from distribution of f(\b|T_i^*>t;\Omega)
  # This distribution can be either normal (using approximation) or t on some user-supplied df.
  
  # Unpack \Omega
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(b.inds, length)); zeta <- O$zeta; sigma <- O$sigma
  b.density <- match.arg(b.density, c('normal', 't'), several.ok = F)
  
  if(b.density == 'normal'){
    # Use optim to draw from f(\b,\Y,T,\Delta;\Omega^{\ell}). Only thing changing per simulation is \Omega
    # Posterior mode and its variance at \Omega^{\ell}
    b.hat.l <- optim(
      b.current, joint_density, joint_density_ddb,
      Y = long$Y, X = long$X, Z = long$Z, beta = beta, D = D, sigma = sigma,
      family = fit$ModelInfo$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
      SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
      b_inds = b.inds, K = length(fit$ModelInfo$family), method = 'BFGS', hessian = T
    )
    Sigma.hat.l <- solve(b.hat.l$hessian)
    b.hat.l <- b.hat.l$par
    # Draw from N(b.hat.l, Sigma.hat.l)
    b.prop.l <- MASS::mvrnorm(n = 1, mu = b.hat.l, Sigma = Sigma.hat.l)
    # DMNVORM on b draws
    current.dens <- as.double(dmvnrm_arma_fast(t(b.current), b.hat.t, Sigma.t, T))
    prop.dens <- as.double(dmvnrm_arma_fast(t(b.prop.l), b.hat.t, Sigma.t, T))
  }else{
    #' Draw from shifted t distribution at \hat{b}, \hat{\Sigma} for subject|T_i>t
    b.prop.l <- mvtnorm::rmvt(n = 1, sigma = Sigma.t, df = df, delta = b.hat.t)
    current.dens <- as.double(dmvt_arma_fast(t(b.current), b.hat.t, Sigma.t, df = df))
    prop.dens <- as.double(dmvt_arma_fast(b.prop.l, b.hat.t, Sigma.t, df = df))
  }
  diff.dens <- current.dens - prop.dens # Difference in current - proposal log-likelihood.
  
  # Joint data log likelihood on current and proposed values of b
  current.joint.ll <- logLik.b(b.current, long, surv, O, beta.inds, b.inds, fit)
  proposed.joint.ll <- logLik.b(c(b.prop.l), long, surv, O, beta.inds, b.inds, fit)
  diff.joint.ll <-  proposed.joint.ll - current.joint.ll
  
  # Accept/reject scheme
  accept <- 0
  a <- exp(diff.joint.ll - diff.dens)
  if(is.nan(a)){
    message("Error in b.mh!!\nInformation on current values:\n\n")
    cat(paste0('jointll.diff:', diff.joint.ll, '\n',
               'dens.diff:', diff.dens,'\n'))
    cat(paste0('b.prop: ', b.prop.l,'\nb.current: ', b.current,'\n',
               'joint b.prop: ', logLik.b(b.prop.l, long, surv, O, beta.inds, b.inds, fit), '\n',
               'joint b.current: ', logLik.b(b.current, long, surv, O, beta.inds, b.inds, fit), '\n'))
  }
  a <- min(1, a)
  U <- runif(1)
  if(U <= a){
    accept <- 1
    b.current <- b.prop.l
  }
  
  list(b.current = c(b.current), accept = accept)
}
