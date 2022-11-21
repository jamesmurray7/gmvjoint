#' @keywords internal
print.summary.joint <- function(x, printD = F, digits = 3){
  if(is.null(fit$SE)) stop('Need to run EM with post.process = T')
  qz <- qnorm(.975)
  .to3dp <- function(x) round(x, 3) # Show all to 3dp only.
  
  # Log-likelihood
  # Print log-likelihood
  cat(.print.loglik(fit))
  cat('\n')
  
  # Model fit info
  K <- length(fit$ModelInfo$ResponseInfo)
  responses <- lapply(sapply(fit$ModelInfo$ResponseInfo, strsplit, '\\s\\('), el, 1)
  families <- unlist(fit$family)
  # Standard errors and parameter estimates.
  SE <- fit$SE
  D <- fit$co$D
  betas <- fit$co$beta
  sigmas <- unlist(fit$co$sigma)
  gammas <- fit$co$gamma
  zetas <- fit$co$zeta
  
  # Random effects matrix
  if(printD){
    cat(paste0('Random effects variance-covariance matrix: \n'))
    print(.to3dp(D))
    cat('\n')
  }
  
  # Longitudinal parts
  MakeTables <- lapply(1:K, function(k){
    
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
    
    cat(paste0(responses[[k]], ' (', families[k], '): \n'))
    my <- c(beta)
    if(sigma != 0.0) my <- c(my, sigma)
    parameter <- names(my)
    
    rSE <- SE[match(names(my), names(SE))]#SE associated with these coeffs
    
    lb <- my - qz * rSE; ub <- my + qz * rSE
    
    z <- my/rSE
    p <- 2 * (pnorm(abs(z), lower.tail = F))
    
    this.out <- setNames(data.frame(.to3dp(my), .to3dp(rSE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                         c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
    print(this.out)
    cat('\n')
  })
  
  # Survival
  cat('Event-time sub-model: \n')
  # Rename gammas?
  survs <- c(zetas, gammas)
  surv.SE <- SE[match(names(survs), names(SE))]
  new.gammas <- sapply(1:K, function(k) gsub('\\_\\d?\\d', paste0('_', unlist(responses)[k]), names(gammas)[k]))
  names(survs)[grepl('gamma', names(survs))] <- new.gammas
  
  lb <- survs - qz * surv.SE; ub <- survs + qz * surv.SE
  z <- survs/surv.SE
  p <- 2 * (pnorm(abs(z), lower.tail = F))
  
  surv.out <- setNames(data.frame(.to3dp(survs), .to3dp(surv.SE), .to3dp(lb), .to3dp(ub), round(p, 3)),
                       c('Estimate', 'SE', '2.5%', '97.5%', 'p-value'))
  print(surv.out)
  
  # Elapsed times
  cat(sprintf('\nElapsed times as follows (%d iterations):\n', fit$iter))
  print(.to3dp(fit$elapsed.time))
  
  invisible(1+1)
}