#' Calculate log-likelihood and degrees of freedom from a multivariate joint model
#' @author James Murray (\email{j.murray7@@ncl.ac.uk})
#' @seealso \code{\link{logLik.joint}}
#' @importFrom mvtnorm dmvnorm
#' @keywords internal
joint.log.lik <- function(coeffs, dmats, b, surv, sv, l0u, l0i, gamma.rep, beta.inds, b.inds, K, q, family, Sigma){
  Y <- dmats$Y; X <- dmats$X; Z <- dmats$Z
  beta <- coeffs$beta; D <- coeffs$D; sigma <- coeffs$sigma; zeta <- coeffs$zeta
  S <- sv$S; SS <- sv$SS; Delta <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi
  beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1)
  
  # log-likelihood conditional on random effects (i.e. not marginalised over REs).
  ll <- mapply(function(Y, X, Z, b, S, SS, Fu, Fi, Delta, l0i, l0u){
    joint_density(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
                  Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
                  beta_inds = beta.inds2, b_inds = b.inds2, K = K) * -1
  }, Y = Y, X = X, Z = Z, b = b, S = S, SS = SS, Fu = Fu, Fi = Fi, Delta = Delta, l0i = l0i, l0u = l0u)
  ll.cond <- sum(ll) 
  
  # Observed data log-likelihood log f(T, Delta, Y|b;Omega)
  ll2 <- mapply(function(Y, X, Z, b, S, SS, Fu, Fi, Delta, l0i, l0u, Sigma){
    
    log.fti <- logfti(b, S, SS, Fi, Fu, l0i, l0u, Delta, gamma.rep, zeta)
    
    log.fY <- lapply(1:K, function(k){
      f <- family[[k]]
      Xk <- X[[k]]; Zk <- Z[[k]]; Yk <- Y[[k]];
      betak <- beta[beta.inds[[k]]]; bk <- b[b.inds[[k]]]
      eta <- Xk %*% betak + Zk %*% bk
      out <- switch(f,
                    gaussian = dnorm(Yk, eta, sqrt(sigma[[k]]), T),
                    binomial = dbinom(Yk, 1, plogis(eta), T),
                    poisson = dpois(Yk, exp(eta), T),
                    genpois = ll_genpois(eta, sigma[[k]], Yk),
                    Gamma = dgamma(Yk, sigma[[k]], exp(eta)/sigma[[k]], log = T))
      out
    })
    log.fY <- sum(do.call(c, log.fY))
    
    log.fb <- 0.5 * q * log(2*pi) + 0.5 * log(det(Sigma)) + mvtnorm::dmvnorm(b, sigma = D, log = T)
    
    c(log.fti, log.fY, log.fb)
    
  }, Y = Y, X = X, Z = Z, b = b, S = S, SS = SS, Fu = Fu, Fi = Fi, Delta = Delta, l0i = l0i, l0u = l0u, Sigma = Sigma)
  ll.obs <- sum(ll2)
  
  # Calculate AIC and BIC
  N <- sum(colSums(do.call(rbind, lapply(Y, function(y) sapply(y, length)))))
  
  df <- sum(dmats$P) + length(vech(D)) +                                   # Fixed effects + D terms,
    sum(unlist(family) %in% c('gaussian', 'genpois', 'Gamma')) +           # Dispersion terms,
    K + ncol(S[[1]])                                                       # K gammas, ncol(S) zetas.
  
  df.residual <- N - (df + sum(dim(do.call(rbind, b))))                    # DF - num. REs
  
  structure(ll.obs,
            'df' = df, 'df.residual' = df.residual,
            'Conditional loglikelihood' = ll.cond,
            'AIC' = -2 * ll.obs + 2 * df,
            'BIC' = -2 * ll.obs + log(N) * df)
}

#' Log-likelihood for joint model.
#' 
#' @description Calculate joint log-likelihood, degrees of freedom, AIC and BIC of 
#' joint model fit.
#' 
#' @param object a \code{joint} object.
#' @param conditional Logical. Should the conditional or observed data log-likelihood
#'  be returned? See \strong{details}.
#' @param ... additional arguments (none used).
#' 
#' @details 
#' 
#' Calculate the log-likelihood of a joint model of survival and multivariate longitudinal
#' data (i.e. a \code{joint} object). The argument \code{conditional} manages whether
#' or not the log-likelihood \emph{conditional} on the random effects, or simply
#' the observed data log-likelihood is returned (the default, \code{conditional = FALSE}). 
#' 
#' If \code{conditional = TRUE}, then the log-likelihood conditional on the random 
#' effects is returned. That is
#' \deqn{\log f(T_i, \Delta_i, Y_i|b_i;\Omega) = 
#'       \log f(Y_i|b_i; \Omega) + \log f(T_i, \Delta_i|b_i; \Omega) + \log f(b_i|\Omega)}
#'       
#' If \code{conditional = FALSE}, then the observed data log-likelihood is returned i.e.
#' 
#' \deqn{\log\int f(Y_i|b_i; \Omega)f(T_i, \Delta_i|b_i; \Omega)f(b_i|\Omega)db_i.}
#' 
#' Additionally, the degrees of freedom, \eqn{nu} is given by
#' 
#' \deqn{\nu = \sum_{k=1}^KP_k + P_s + \code{length(vech(D))} + P_{\sigma_k},}
#' 
#' where \eqn{P_k} denotes the number of coefficients estimated for the \eqn{k}th response,
#' and \eqn{P_{\sigma_k}} the number of dispersion parameters estimated. \eqn{P_s} denotes
#' the number of survival coefficients, i.e. the length of \code{c(zeta, gamma)}. Finally,
#' all covariance parameters are captured in \code{length(vech(D))}. 
#' 
#' With the degrees of freedom, we can additionally compute AIC and BIC, which are defined
#' in no special way; and are calculated using the observed data log-likelihood.
#' 
#' @seealso \code{\link{joint.log.lik}}
#' @author James Murray (\email{j.murray7@@ncl.ac.uk})
#' 
#' @references 
#'  
#' Henderson R, Diggle P, Dobson A. Joint modelling of longitudinal measurements and event time
#' data. \emph{Biostatistics} 2000; \strong{1(4)}; 465-480.
#' 
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data measured with error.
#' \emph{Biometrics} 1997; \strong{53(1)}; 330-339.
#' 
#' @keywords methods
#' @export
logLik.joint <- function(object, conditional = FALSE, ...){
  if(!inherits(object, 'joint')) stop("Only usable with object of class 'joint'.")
  if(is.null(object$logLik)) stop("Rerun with post.process = TRUE.")
  ll <- object$logLik
  class(ll) <- 'logLik'
  ll
}

##' Extract AIC from a joint model fit.
##'
##' @param fit A fitted \code{joint} object,
##' @param scale See \code{\link[stats]{extractAIC}}; not used.
##' @param k Numeric specifying the "weight" of degrees of freedom (default \code{k=2}).
##' @param conditional Should AIC of conditional or observed log-likelihood be used? Defaults
##' to \code{conditional = FALSE}.
##' @param ... additional arguments (none used).
##'
##' @method extractAIC joint
##' @export
extractAIC.joint <- function(fit, scale, k = 2, conditional = FALSE, ...){
  x <- fit
  if(!inherits(x, 'joint')) stop("Only usable with object of class 'joint'.")
  if(is.null(x$logLik)) stop("Rerun with post.process = TRUE.")
  L <- x$logLik
  df <- attr(L, 'df')
  if(conditional){
    ll <- c(attr(L, 'Conditional loglikelihood'))
  }else{
    ll <- c(L)
  }
  return(setNames(c(df, c(-2*ll + k * df)), c('df', 'AIC')))
}




