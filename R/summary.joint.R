#' Summary of an \code{joint} object.
#' 
#' @description Generate summary of a fitted multivariate joint model.
#' 
#' @author James Murray \email{j.murray7@@ncl.ac.uk}
#' @method summary joint
#' @export
#' @seealso \code{\link{joint}} and \code{\link{joint.object}}
#' 
#' @param object a joint model fit by the \code{joint} function.
#' @param ... additional arguments (none used).
#' 
#' @returns Object of class \code{summary.joint}.
#'  
#' @examples 
#' \donttest{
#' data(PBC)
#' long.formula <- list(
#'   platelets ~ time * drug + (1 + time|id),
#'   albumin ~ time * drug + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ sex + drug
#'
#' PBC <- na.omit(PBC[,c('id', 'survtime', 'status', 'sex', 
#'                       'drug', 'platelets', 'albumin', 'time')])
#' fit <- joint(long.formula, surv.formula, PBC, family = list('genpois', 'gaussian'),
#'             control = list(verbose = TRUE))
#' summary(fit)
#' }
summary.joint <- function(object, ...){
  if(!inherits(object, 'joint')) stop("Only usable with object of class 'joint'.")
  if(is.null(object$SE)) stop('Rerun with control(post.process = TRUE).')
  qz <- qnorm(.975)
  
  # Extract things to ModelInfo
  M <- object$ModelInfo
  K <- M$K                                         # Number of responses
  responses <- M$Resps                             # Response names only
  families <- unlist(M$family)                     # Families only
  nobs <- M$nobs; n <- M$n; nev <- M$nev           # Data stuff
  long.formulas <- M$long.formulas                 # Sub-models
  surv.formula <- M$surv.formula 
  disp.formulas <- M$disp.formulas                 # (+ dispersion)
  inds.beta <- M$inds$R$beta; inds.b <- M$inds$R$b # indices
  
  # Parameter estimates and standard errors.
  SE <- object$SE
  coeffs <- object$coeffs
  D <- coeffs$D
  betas <- coeffs$beta
  sigmas <- coeffs$sigma
  gammas <- coeffs$gamma
  zetas <- coeffs$zeta
  
  # Longitudinal parts
  Longits <- setNames(lapply(1:K, function(k){
    
    beta <- betas[inds.beta[[k]]]
    beta.SE <- SE[match(names(beta), names(SE))]
    beta.lb <- beta - qz * beta.SE; beta.ub <- beta + qz * beta.SE
    
    if(families[k] %in% c("gaussian", "genpois", "Gamma", "negbin")){
      sigma <- sigmas[[k]]
      sigma.SE <- SE[match(names(sigma), names(SE))]
      sigma.lb <- sigma - qz * sigma.SE; sigma.ub <- sigma + qz * sigma.SE
      if(families[k] %in% c("Gamma", "negbin")){
        sigma <- exp(sigma)
        sigma.SE <- sigma.SE * sigma
        sigma.lb <- exp(sigma.lb); sigma.ub <- exp(sigma.ub)
      }
    }else{
      sigma <- sigma.SE <- sigma.lb <- sigma.ub <- NULL
    }
    
    beta <- c(beta, sigma); rSE <- c(beta.SE, sigma.SE); lb <- c(beta.lb, sigma.lb); ub <- c(beta.ub, sigma.ub)
    parameter <- names(beta)
    
    z <- beta/rSE
    p <- 2 * (pnorm(abs(z), lower.tail = F))
    
    this.out <- setNames(data.frame(beta, rSE, z, p, lb, ub),
                         c('Estimate', 'SE', "Z", "p-value", '2.5%', '97.5%'))
    this.out
  }), unlist(M$ResponseInfo))
  
  # Survival parts
  survs <- c(zetas, gammas)
  surv.SE <- SE[match(names(survs), names(SE))]
  new.gammas <- sapply(1:K, function(k) gsub('\\_\\d?\\d', paste0('_', unlist(responses)[k]), names(gammas)[k]))
  names(survs)[grepl('gamma', names(survs))] <- new.gammas
  
  lb <- survs - qz * surv.SE; ub <- survs + qz * surv.SE
  z <- survs/surv.SE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  Survs <- setNames(data.frame(survs, surv.SE, z, p, lb, ub),
                       c('Estimate', 'SE', 'Z', 'p-value', '2.5%', '97.5%'))
  
  # Other items to return
  et <- object$elapsed.time; iter <- unname(as.integer(et['iterations'])); et <- et[-which(names(et) == 'iterations')]
  haz <- object$hazard
  ll <- object$logLik
  
  out <- list(
    # Model stuff
    responses = responses,
    families = families,
    long.formulas = long.formulas, surv.formula = surv.formula,
    nobs = nobs, n = n, nev = nev,
    # Coefficients
    Longits = Longits, Survs = Survs, haz = haz, SE = SE, D = D,
    et = et, iter = iter, ll = ll
  )
  class(out) <- 'summary.joint'
  out
}

#' @keywords internal
#' @method print summary.joint
#' @export
print.summary.joint <- function(x, digits = 3, printD = FALSE, ...){
  if(!inherits(x, 'summary.joint')) stop("Only usable with object of class 'summary.joint'.")
  .round <- function(X) round(X, digits) # Show all to <digits> only.
  K <- length(x$families)
  
  # Print data information
  cat(sprintf("\nNumber of subjects: %d\n", x$n))
  cat(sprintf("Number of events: %d (%.2f%%)\n", x$nev, 100 * x$nev/x$n))
  
  # Print log-likelihood
  cat("\nModel fit statistics ----\n")
  ll <- x$ll
  print(c("log.Lik" = c(ll),
          "df" = c(attr(ll, 'df')),
          "AIC" = c(attr(ll, 'AIC')),
          "BIC" = c(attr(ll, 'BIC'))))
  
  # Print variance-covariance matrix on REs
  if(printD){
    cat("\nRandom effects variance-covariance ----\n")
    print(.round(x$D))
  }
  
  if(K > 1) cat("\nLongitudinal processes ----") else cat("\nLongitudinal Process ----")
  
  # Print K coefficients for longitudinal estimates.
  lapply(1:K, function(k){
    cat(sprintf("\n%s (%s)\nCall:\n", x$responses[[k]], x$families[k]))
    cat(deparse(x$long.formulas[[k]]), '\n')
    L <- x$Longits[[k]]
    print(cbind(L[,1], apply(L[,-1], 2, .round)))
  })
  
  # Survival coefficients
  cat(sprintf('\nEvent-time sub-model: ----\nCall: %s\n', deparse(x$surv.formula)))
  print(cbind(x$Survs[,1], apply(x$Survs[,-1],2,.round)))
  
  
  # Computation summary
  cat('\nComputation summary ----\n')
  cat(sprintf("Number of EM iterations: %d,\nTime spent in EM algorithm: %.2fs\nTotal computation time: %.2fs.\n", 
              x$iter, x$et['EM time'], x$et['Total Computation time']))
  
  cat("\n")
  invisible(x)
}


