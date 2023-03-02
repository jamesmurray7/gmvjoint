# Half-vectorisation of matrix x
#' @keywords internal
vech <- function(x) x[lower.tri(x, T)]

# Parsing input formula
#' @keywords internal
#' @importFrom glmmTMB splitForm
parseFormula <- function(formula){ 
  split <- glmmTMB::splitForm(formula, allowFixedOnly = F)
  fixed <- split$fixedFormula
  random <- el(split$reTrmFormulas)
  
  # Parse fixed effects
  response <- as.character(fixed)[2]
  fixed <- as.character(fixed)[3]
  if(grepl('splines\\:\\:|ns\\(|bs\\(', fixed)){
    attr(fixed, 'special') <- 'spline'
    if(grepl('ns\\(', fixed)) attr(fixed, 'spline type') <- 'natural'
    else if(grepl('bs\\(', fixed)) attr(fixed, 'spline type') <- 'basis'
    else stop('Unknown spline type')
  }else{
    attr(fixed, 'special') <- 'none'
  }
  
  # Parse random effects
  random.by <- as.character(random)[3]
  random <- as.character(random)[2]
  if(grepl('splines\\:\\:|ns\\(|bs\\(', random)){
    attr(random, 'special') <- 'spline'
    if(grepl('ns\\(', random)) attr(random, 'spline type') <- 'natural'
    else if(grepl('bs\\(', random)) attr(random, 'spline type') <- 'basis'
    else stop('Unknown spline type')
  }else{
    attr(random, 'special') <- 'none'
  }
  
  return(list(
    response = response,
    fixed = fixed,
    random = random,
    random.by = random.by
  ))
}

#' @keywords internal
converge.check <- function(params.old, params.new, criteria, iter, Omega, verbose){
  
  type <- criteria$type
  # Absolute difference
  diffs.abs <- abs(params.new - params.old)
  # Relative difference
  diffs.rel <- diffs.abs/(abs(params.old) + criteria$tol.den)
  # SAS convergence criterion
  sas.crit <- abs(params.old) >= criteria$threshold
  sas.abs <- diffs.abs < criteria$tol.abs
  sas.rel <- diffs.rel < criteria$tol.rel
  sas.conv <- all(sas.abs[!sas.crit]) & all(sas.rel[sas.crit])
  
  # Check convergence based on user-supplied criterion
  if(type == "abs"){
    converged <- max(diffs.abs) < criteria$tol.abs
  }else if(type == "rel"){
    converged <- max(diffs.rel) < criteria$tol.rel
  }else if(type == "either"){
    converged <- (max(diffs.abs) < criteria$tol.abs) | (max(diffs.rel) < criteria$tol.rel)
  }else if(type == "sas"){
    converged <- sas.conv
  }
  
  if(verbose){
      cat("\n")
      cat(sprintf("Iteration %d:\n", iter))
      cat("vech(D):", round(vech(Omega$D), 4), "\n")
      cat("beta:", round(Omega$beta, 4), "\n")
      if(any(unlist(Omega$sigma) != 0)) cat("sigma:", round(unlist(Omega$sigma)[unlist(Omega$sigma) != 0], 4), "\n")
      cat("gamma:", round(Omega$gamma, 4), "\n")
      cat("zeta:", round(Omega$zeta, 4), "\n")
      cat("\n")
      cat("Maximum absolute difference:", round(max(diffs.abs), 4), "for",
          names(params.new)[which.max(diffs.abs)], "\n")
      cat("Maximum relative difference:", round(max(diffs.rel), 4), "for",
                 names(params.new)[which.max(diffs.rel)], "\n")
      if(converged) cat("Converged! (Criteria:", type, ").\n\n")
  }
  
  list(converged = converged,
       diffs.abs = diffs.abs, diffs.rel = diffs.rel)
}

# Don't think this is ever used  -- remove?
# Create appropriately-dimensioned matrix of random effects.
#' @keywords internal
bind.bs<- function(bsplit){
  qmax <- max(sapply(bsplit, length)) # Maximum number of REs across all longitudinal responses.
  # Pad with zeros until each row is of length qmax
  step <- lapply(bsplit, function(b){
    l <- length(b)
    if(l<qmax) b <- c(b, rep(0, (qmax-l)))
    b
  })
  step <- do.call(rbind, step); colnames(step) <- NULL
  as.matrix(step)
}

# Not used, keeping anyway.
# Obtain hessian from a score vector using some differencing method.
#' @keywords internal
numDiff <- function(x, f, ..., method = 'central', heps = 1e-4){
  method <- match.arg(method, c('central', 'forward', 'Richardson'))
  n <- length(x)
  out <- matrix(0, nrow = n, ncol = n)
  hepsmat <- diag(pmax(abs(x), 1) * heps, nrow = n)
  if(method == "central"){
    for(i in 1:n){
      hi <- hepsmat[,i]
      fdiff <- c(f(x + hi, ...) - f(x - hi, ...))
      out[, i] <- fdiff/(2 * hi[i])
    }
  }else if(method == "Richardson"){
    for(i in 1:n){
      hi <- hepsmat[,i]
      fdiff <- c(f(x - 2 * hi, ...) - 8 * f(x - hi, ...) + 8 * f(x + hi, ...) - f(x + 2 * hi, ...))
      out[,i] <- fdiff/(12*hi[i])
    }
  }else{
    f0 <- f(x, ...)
    for(i in 1:n){
      hi <- hepsmat[,i]
      fdiff <- c(f(x + hi, ...) - f0)
      out[,i] <- fdiff/hi[i]
    }
  }
  (out + t(out)) * .5
}
