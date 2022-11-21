#' Calculate joint log-likelihood, degrees of freedom, AIC and BIC of joint model fit.
#' 
#' @references 
#' 
#' Henderson R, Diggle P, Dobson A. Joint modelling of longitudinal measurements and event time data.
#' \emph{Biostatistics} 2000; \strong{1(4)}; 465-480.
#' 
#' @keywords internal
joint.log.lik <- function(coeffs, dmats, b, surv, sv, l0u, l0i, gamma.rep, beta.inds, b.inds, K, q, family){
  Y <- dmats$Y; X <- dmats$X; Z <- dmats$Z
  beta <- coeffs$beta; D <- coeffs$D; sigma <- coeffs$sigma; zeta <- coeffs$zeta
  S <- sv$S; SS <- sv$SS; Delta <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi
  beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) # Indexed for C++ use.
  # f(Y|b; \Omega) + f(T|\Omega) + f(\b|D) ---------------
  # (This is identical to joint density) with no contribution from \b in logf(T|...).
  ll1 <- mapply(function(Y, X, Z, b, S, SS, Fu, Fi, Delta, l0i, l0u){
    joint_density(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
                  Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = rep(0, length(gamma.rep)), zeta = zeta,
                  beta_inds = beta.inds2, b_inds = b.inds2, K = K) * -1
  }, Y = Y, X = X, Z = Z, b = b, S = S, SS = SS, Fu = Fu, Fi = Fi, Delta = Delta, l0i = l0i, l0u = l0u)
  out <- sum(ll1)
  
  # Calculate AIC and BIC
  N <- sum(colSums(do.call(rbind, lapply(Y, function(y) sapply(y, length)))))
  
  df <- sum(dmats$P) + length(vech(D)) +                                   # Fixed effects + D terms,
    sum(unlist(family) %in% c('gaussian', 'genpois', 'Gamma')) +       # Dispersion terms,
    K + ncol(S[[1]])                                                   # K gammas, ncol(S) zetas.
  
  df.residual <- N - (df + sum(dim(do.call(rbind, b))))                    # DF - num. REs
  
  structure(out,
            'df' = df, 'df.residual' = df.residual,
            'AIC' = -2 * out + 2 * df,
            'BIC' = -2 * out + log(N) * df)
}


#' @keywords internal
logLik.joint <- function(fit){
  ll <- fit$logLik
  df <- attr(ll, 'df')
  aic <- attr(ll, 'AIC')
  cat(sprintf("logLik: %.3f (df: %d), AIC: %.3f\n", ll, df, aic))
}
