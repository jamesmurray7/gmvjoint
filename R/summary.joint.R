#' Summary of an \code{joint} object.
#' 
#' @description Generate summary of a fitted multivariate joint model.
#' 
#' @author James Murray \email{j.murray7@@ncl.ac.uk}
#' @keywords internal
#' @method summary joint
#' @seealso \code{\link{joint}} and \code{\link{joint.object}}
#' 
#' @param x a joint model fit by the \code{joint} function.
#' 
#' @returns Object of class \code{summary.joint}.
#' 
#' @export
#' 
#' @examples 
#' data(PBC)
#' long.formula <- list(
#'   platelets ~ time * drug + (1 + time|id),
#'   albumin ~ time * drug + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ sex + drug
#'
#' PBC <- na.omit(PBC[,c('id', 'survtime', 'status', 'sex', 
#' 'drug', 'platelets', 'albumin', 'time')])
#' fit <- joint(long.formula, surv.formula, PBC, family = list('genpois', 'gaussian'),
#'             control = list(verbose = TRUE))
#' summary(fit)
summary.joint <- function(x, ...){
  if(!inherits(x, 'joint')) stop("Only usable with object of class 'joint'.")
  if(is.null(x$SE)) stop('Rerun with post.process = TRUE.')
  qz <- qnorm(.975)
  
  # Extract things to ModelInfo
  M <- x$ModelInfo
  K <- length(M$ResponseInfo)                      # Number of responses
  responses <- lapply(sapply(M$ResponseInfo,       # Response names only
                             strsplit, '\\s\\('), el, 1)
  families <- unlist(M$family)                     # Families only
  nobs <- M$nobs; n <- M$n; nev <- M$nev           # Data stuff
  long.formulas <- M$long.formulas                 # Sub-models
  surv.formula <- M$surv.formulas
  inds.beta <- M$inds$beta; inds.b <- M$inds$b     # indices
  
  # Parameter estimates and standard errors.
  SE <- x$SE
  coeffs <- x$coeffs
  D <- coeffs$D
  betas <- coeffs$beta
  sigmas <- unlist(coeffs$sigma)
  gammas <- coeffs$gamma
  zetas <- coeffs$zeta
  
  # Longitudinal parts
  Longits <- setNames(lapply(1:K, function(k){
    
    beta <- betas[grepl(responses[[k]], names(betas))]
    sigma <- setNames(sigmas[k], 
                      if(families[k] == "genpois")
                        paste0("phi_", k)
                      else if(families[k] == 'Gamma')
                        paste0("shape_", k)
                      else if(families[k] == 'gaussian')
                        paste0('sigma^2_', k)
                      else
                        'NO DISP')
    
    betak <- c(beta)
    if(sigma != 0.0) betak <- c(betak, sigma)
    parameter <- names(betak)
    
    rSE <- SE[match(names(betak), names(SE))]#SE associated with these coeffs
    
    lb <- betak - qz * rSE; ub <- betak + qz * rSE
    
    z <- betak/rSE
    p <- 2 * (pnorm(abs(z), lower.tail = F))
    
    this.out <- setNames(data.frame(betak, rSE, z, p, lb, ub),
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
  et <- x$elapsed.time; iter <- unname(as.integer(et['iterations'])); et <- et[-which(names(et) == 'iterations')]
  haz <- x$hazard
  ll <- x$logLik
  
  
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

#' @method print summary.joint
#' @keywords internal
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
  invisible(1+1)
}


