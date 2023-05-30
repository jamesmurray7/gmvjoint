#' @importFrom mvtnorm dmvnorm
#' @keywords internal
joint.log.lik <- function(coeffs, dmats, surv, sv, family, b, l0i, l0u, inds, Sigma){
  beta <- coeffs$beta; D <- coeffs$D; sigma <- coeffs$sigma; zeta <- coeffs$zeta; gamma <- coeffs$gamma
  gamma.rep <- rep(gamma, sapply(inds$R$b, length))
  
  # log-likelihood conditional on random effects (i.e. not marginalised over REs).
  ll <- mapply(function(Y, X, Z, W, b, S, SS, Fu, Fi, Delta, l0i, l0u){
    joint_density(b = b, Y = Y, X = X, Z = Z, W = W, beta = beta, D = D, sigma = sigma, family = family,
                  Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
                  beta_inds = inds$Cpp$beta, b_inds = inds$Cpp$b, K = dmats$K) * -1
  }, Y = dmats$Y, X = dmats$X, Z = dmats$Z, W = dmats$W,
     b = b, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, 
     Delta = surv$Delta, l0i = l0i, l0u = l0u)
  ll.cond <- sum(ll) 
  
  # Observed data log-likelihood log f(T, Delta, Y|b;Omega)
  ll.obs <- ll + 0.5 * log(2*pi) * sv$q + sapply(Sigma, function(x) 0.5 * log(det(x)))
  ll.obs <- sum(ll.obs)
  
  # Calculate AIC and BIC
  N <- sum(dmats$m)
  
  df <- sum(dmats$P) + length(vech(D)) +   # Fixed effects + D terms,
        sum(unlist(sigma) != 0L) +         # Dispersion
        dmats$K + ncol(sv$S[[1]])          # K gammas, ncol(S) zetas.
  
  # Residual df: total nobs - df - number of random effects
  df.residual <- N - (df + sum(dim(do.call(rbind, b))))                    
  
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
#' Additionally, the degrees of freedom, \eqn{\nu} is given by
#' 
#' \deqn{\nu = \code{length(vech(D))} + \sum_{k=1}^K\{P_k + P_{\sigma_k}\} + P_s,}
#' 
#' where \eqn{P_k} denotes the number of coefficients estimated for the \eqn{k}th response,
#' and \eqn{P_{\sigma_k}} the number of dispersion parameters estimated. \eqn{P_s} denotes
#' the number of survival coefficients, i.e. the length of \code{c(zeta, gamma)}. Finally,
#' all covariance parameters are captured in \code{length(vech(D))}. 
#' 
#' With the degrees of freedom, we can additionally compute AIC and BIC, which are defined
#' in no special way; and are calculated using the observed data log-likelihood.
#' 
#' @seealso \code{\link{extractAIC.joint}} and \code{\link{anova.joint}}
#' @author James Murray (\email{j.murray7@@ncl.ac.uk})
#' 
#' @returns Returns an object of class \code{logLik}, a number which is the log-likelihood
#' of the fitted model \code{object}. This has multiple attributes: \code{df} which is the 
#' degrees of freedom, \code{df.residual}; the number of residual degrees of freedom;
#' \code{AIC} and {BIC} which are the Akaike or Bayes information criterion evaluated at either
#' the conditional or observed log-likelihood (as requested by argument \code{conditional}).
#' 
#' @examples 
#' # Bivariate simulated data (2x Gaussian)
#' data <- simData(n = 100,
#'    D = diag(c(.25, .04, .2, .02)),
#'    gamma = c(0.4, -0.2), theta = c(-2, .2))$data
#' fit <- joint(list(
#'     Y.1 ~ time + cont + bin + (1 + time|id),
#'     Y.2 ~ time + cont + bin + (1 + time|id)
#'   ), Surv(survtime, status) ~ cont + bin, 
#'   data = data, 
#'   family = list('gaussian', 'gaussian'))
#' 
#' logLik(fit)
#' 
#' @references 
#'  
#' Henderson R, Diggle P, Dobson A. Joint modelling of longitudinal measurements and event time
#' data. \emph{Biostatistics} 2000; \strong{1(4)}; 465-480.
#' 
#' Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data measured with error.
#' \emph{Biometrics} 1997; \strong{53(1)}; 330-339.
#' 
#' @method logLik joint
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
##' 
##' @returns A numeric vector of length 2, with first and second element giving \describe{
##' \item{\code{df}}{The degrees of freedom for the fitted model.}
##' \item{\code{AIC}}{The Akaike Information Criterion for the fitted model.}
##' }
##' 
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


#' Anova for joint models
#' 
#' @description Perform a likelihood ratio test between two (nested) \code{joint} models. 
#'
#' @param object a joint model fit by the \code{joint} function. This should be \strong{nested}
#' in \code{object2}.
#' @param object2 a joint model fit by the \code{joint} function. This should be more complex
#' than \code{object} whilst sharing the same survival sub-model.
#' @param ... additional arguments (none used).
#'
#' @return A list of class \code{anova.joint} with elements \describe{
#' 
#'   \item{\code{mod0}}{the name of \code{object}.}
#'   \item{\code{l0}}{the log-likelihood of the nested model, i.e. fit under the null.}
#'   \item{\code{AIC0}}{AIC for \code{object}.}
#'   \item{\code{BIC0}}{BIC for \code{object}.}
#'   \item{\code{mod1}}{the name of \code{object2}.}
#'   \item{\code{l1}}{the log-likelihood under the alternative hypothesis.}
#'   \item{\code{AIC1}}{AIC for \code{object2}.}
#'   \item{\code{BIC1}}{BIC for \code{object2}.}
#'   \item{\code{LRT}}{likelihood ratio test statistic.}
#'   \item{\code{p}}{the p-value of \code{LRT}.}
#' 
#' }
#' @export
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk})
#' @seealso \code{\link{joint}} and \code{\link{logLik.joint}}.
#' @method anova joint
#'
#' @examples
#' \donttest{
#' rm(list=ls())
#' data(PBC)
#' # Compare quadratic vs linear time specification for log(serum bilirubin) -----
#' PBC$serBilir <- log(PBC$serBilir)
#' long.formulas1 <- list(serBilir ~ drug * time + (1 + time|id))
#' long.formulas2 <- list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' # Fit the two competing models (fit is nested in fit2) ------------------------
#' fit <- joint(long.formulas1, surv.formula, PBC, family, 
#'              control = list(verbose = FALSE))
#' fit2 <- joint(long.formulas2, surv.formula, PBC, family, control = list(verbose = FALSE))
#' anova(fit, fit2)
#' # Quadratic terms improve fit significantly. 
#' }
anova.joint <- function(object, object2, ...){
  if(missing(object) || missing(object2))
    stop("both 'object' and 'object2' should be provided (single-object anova not yet implemented).")
  if(!inherits(object,  'joint')) stop("Only usable with object of class 'joint'.")
  if(!inherits(object2, 'joint')) stop("Only usable with object of class 'joint'.")
  

  # Ensure all constituent families are the same (?)
  K0 <- length(object$ModelInfo$family); K1 <- length(object2$ModelInfo$family)
  if(K0 != K1)
    warning("Comparison between two models of different dimension.")
  if(K0 == K1){
    F0 <- unlist(object$ModelInfo$family); F1 <- unlist(object2$ModelInfo$family)
    check <- all(F0 == F1)
    if(!check) stop("Family mismatch in model specifications for '", deparse(substitute(object)), "' and '",
                    deparse(substitute(object)), "'.")
  }
  
  # Check models were fit to same data
  if(object$ModelInfo$nobs!=object2$ModelInfo$nobs) stop("Models were not fit to same data.")
  
  # Extract lls
  l0 <- logLik(object); l1 <- logLik(object2)
  # Extract dfs and ensure object2 more complex.
  df0 <- attr(l0, 'df'); df1 <- attr(l1, 'df')
  xdf <- df1 - df0
  if(xdf < 0)
    stop("'", deparse(substitute(object)), "' should be nested in '", deparse(substitute(object2)), "' (i.e. '",
         deparse(substitute(object2)), "' more complex).")
  
  LRT <- -2 * (c(l0) - c(l1))
  p <- pchisq(LRT, xdf, lower.tail = F)
  
  out <- list(
    mod0 = deparse(substitute(object)),
    l0   = c(l0), AIC0 = AIC(object), BIC0 = attr(l0, 'BIC'), df0 = df0,
    mod1 = deparse(substitute(object2)),
    l1   = c(l1), AIC1 = AIC(object2), BIC1 = attr(l1, 'BIC'), df1 = df1,
    LRT = LRT, p = p
  )
  class(out) <- 'anova.joint'
  out
}

#' @method print anova.joint
#' @keywords internal
#' @export
print.anova.joint <- function(x, ...){
  if(!inherits(x, 'anova.joint')) stop("Only usable with objects of class 'anova.joint'.")
  p <- x$p
  if(p < 1e-3) p <- '< 0.001' else p <- as.character(round(p, 3))
  df <- setNames(data.frame(logLik = c(x$l0, x$l1), df = c(x$df0, x$df1),
                 AIC = c(x$AIC0, x$AIC1), BIC = c(x$BIC0, x$BIC1),
                 lrt = c("", as.character(round(x$LRT, 3))),
                 p = c("", p)),
                 c("logLik", "df", "AIC", "BIC", "Likelihood ratio test", "p-value"))
  row.names(df) <- c(x$mod0, x$mod1)
  cat("\n")
  print(df)
  cat("\n")
}

