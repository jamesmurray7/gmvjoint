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

#' Take and return difference between two vectors according to some criterion.
#' @keywords internal
difference <- function(params.old, params.new, type){
  if(type == 'absolute'){
    rtn <- abs(params.new - params.old)
  }else if(type == 'relative'){
    rtn <- max(
      (params.new - params.old)/(abs(params.old) + 1e-3)
    )
  }else{
    rtn <- NA
  }
  rtn
}

# Don't think this is ever used (?)
#' create appropriately-dimensioned matrix of random effects.
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
