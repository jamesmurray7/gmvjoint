#' @keywords internal
#' @importFrom MASS mvrnorm
#' @importFrom pracma nearest_spd
#' @importFrom stats vcov
Omega.draw <- function(x){
  Omega.Var <- vcov(x)
  # Mean vector
  co <- x$coeffs
  sigma <- unlist(co$sig)
  sigma <- sigma[sigma != 0L]
  Omega.mean <- setNames(
    c(
      c(vech(co$D)),
      co$beta, sigma,
      co$gamma, co$zeta
    ), names(x$SE))
  
  # Draw from N(Omega.mean, Omega.Var)
  non.invertible <- TRUE            # Temporary fix --> ensure we get a variance covariance matrix which is invertible
  attempts <- 0
  while(non.invertible){
    draw <- MASS::mvrnorm(n = 1, mu = Omega.mean, Sigma = Omega.Var)
    # Re-construct Omega at this draw.
    D <- vech2mat(x = draw[grepl("^D\\[", names(draw))], x$ModelInfo$Pcounts$q)
    # Check this is pos-definite and transform if not.
    if(is.not.SPD(D)) D <- pracma::nearest_spd(D)
    non.invertible <- inherits(try(solve(D), silent = T), 'try-error')
    attempts <- attempts + 1
  }
  
  
  beta <- draw[match(names(co$beta), names(draw))]
  gamma <- draw[match(names(co$gamma), names(draw))]
  zeta <- draw[match(names(co$zeta), names(draw))]
  # Sort out sigma --> MUST be a list.
  .sigma <- lapply(co$sigma, function(x){
    if(length(x) == 1L & x == 0L)
      return(0)
    else
      return(draw[match(names(x), names(draw))])
  })
  
  list(
    D = D, beta = beta, sigma = .sigma, gamma = gamma, zeta = zeta
  )
}

# Not used 24/10/23 -->
# Complete data log-likelihood w.r.t random effects b.
#' @keywords internal
logLik_b <- function(b, long, surv, O, beta.inds, b.inds, fit){
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(b.inds, length)); zeta <- O$zeta; sigma <- O$sigma
  
  neg.ll <- joint_density(b = b, Y = long$Y, X = long$X, Z = long$Z, W = long$W, 
                          beta = beta, D = D, sigma = sigma,
                          family = fit$ModelInfo$family, Delta = surv$Delta, S = surv$S, Fi = surv$Fi, l0i = surv$l0i,
                          SS = surv$SS, Fu = surv$Fu, haz = surv$l0u, gamma_rep = gamma, zeta = zeta, beta_inds = beta.inds,
                          b_inds = b.inds, K = fit$ModelInfo$K)
  
  -neg.ll
}

# Metropolis-Hastings scheme to draw random effects from candidate density.
#' @keywords internal
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
b.mh <- function(b.current, b.hat.t, Sigma.t, long, surv, O, fit, df){
  # Metropolis-Hasting scheme to draw from distribution of f(\b|T_i^*>t;\Omega)
  # This distribution can be either normal (using approximation) or t on some user-supplied df.
  
  # Unpack \Omega^{ell}
  beta <- O$beta; D <- O$D; gamma <- rep(O$gamma, sapply(fit$ModelInfo$inds$Cpp$b, length)); zeta <- O$zeta; sigma <- O$sigma
  #' Draw from shifted t distribution at \hat{b}, \hat{\Sigma} for subject|T_i>t
  b.prop.l <- c(mvtnorm::rmvt(n = 1, sigma = Sigma.t, df = df, delta = b.hat.t))
  
  # Joint data log likelihood on current and proposed values of b
  logf.current <- -1 * joint_density(b.current, long$Y, long$X, long$Z, long$W, beta,
                                     D, sigma, fit$ModelInfo$family, 0L, surv$S, surv$Fi,
                                     0., surv$SS, surv$Fu, surv$l0u, gamma, zeta, 
                                     fit$ModelInfo$inds$Cpp$beta, fit$ModelInfo$inds$Cpp$b, fit$ModelInfo$K)
  logf.proposal <- -1 * joint_density(b.prop.l, long$Y, long$X, long$Z, long$W, beta,
                                      D, sigma, fit$ModelInfo$family, 0L, surv$S, surv$Fi,
                                      0., surv$SS, surv$Fu, surv$l0u, gamma, zeta, 
                                      fit$ModelInfo$inds$Cpp$beta, fit$ModelInfo$inds$Cpp$b, fit$ModelInfo$K)
  log.firstbit <- logf.proposal - logf.current
  
  log.dens.current <- as.double(dmvt_fast(b.current, b.hat.t, Sigma.t, df))
  log.dens.proposal <- as.double(dmvt_fast(b.prop.l, b.hat.t, Sigma.t, df))
  log.secondbit <- log.dens.current - log.dens.proposal
  
  # Accept/reject scheme
  accept <- 0
  a <- min(exp(log.firstbit - log.secondbit), 1.)
  if(is.nan(a)){
    message("\n\nError in b.mh\nInformation on current values:")
    cat("\n--------\n")
    cat("Current value of b:", round(b.current, 3), '\n')
    cat("Value log f(b current|D):", round(log.dens.current, 3), "\n")
    cat("Current log joint density:", round(logf.current))
    cat("\n--------\n")
    cat("Proposal value of b:", round(b.prop.l, 3), '\n')
    cat("Value log f(b proposal|D): ", round(log.dens.proposal, 3), "\n")
    cat("Proposal log joint density:", round(logf.proposal))
    cat("\n--------\n")
    tt <- try(solve(O$D), silent = TRUE)
    if(inherits(tt, 'try-error')) cat("Generated covariance matrix D non-invertible --> likely issue.\n")
  }
  U <- runif(1)
  if(U <= a){
    accept <- 1
    b.current <- b.prop.l
  }
  
  list(b.current = c(b.current), accept = accept)
}
