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

# Take and return difference between two vectors according to some criterion.
difference <- function(params.old, params.new, type){
  if(type == 'absolute'){
    rtn <- abs(params.new - params.old)
  }else if(type == 'relative'){
    rtn <- abs(params.new - params.old)/(abs(params.old) + 1e-3)
  }else{
    rtn <- NA
  }
  rtn
}

#' @keywords internal
converge.check <- function(params.old, params.new, criteria, iter, Omega, verbose){
  
  type <- criteria$type
  # Absolute difference
  diffs.abs <- abs(params.new - params.old)
  # Relative difference
  diffs.rel <- abs(params.new - params.old)/(abs(params.old) + criteria$tol.den)
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
      cat(paste0("Maximum absolute difference: ", round(max(diffs.abs), 4), " for ",
          names(params.new)[which.max(diffs.abs)], "\n"))
      cat(paste0("Maximum relative difference: ", round(max(diffs.rel), 4), " for ",
                 names(params.new)[which.max(diffs.rel)], "\n"))
      if(converged) cat(paste0("Converged!\n\n"))
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

# Obtain hessian from a score vector using central difference
# NB not used anywhere -- to remove.
#' @keywords internal
cendiff <- function(x, f, ..., eps = .Machine$double.eps^(1/4)){
  n <- length(x)
  out <- matrix(0, nrow = n, ncol = n)
  xi <- pmax(abs(x), 1) * eps
  for(i in 1:n){
    a <- b <- x
    a[i] <- x[i] + xi[i]
    b[i] <- x[i] - xi[i]
    fdiff <- c(f(a, ...) - f(b, ...))
    xdiff <- a[i] - b[i]
    out[, i] <- fdiff/xdiff
  }
  (out + t(out)) * .5
}
