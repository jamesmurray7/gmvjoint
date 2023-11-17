# For bootstrapping purposes makes sense to take in a set of parameter 
# initial conditions, or at least offer this anyway (e.g. known values).
# The point of this is to start EM "close" to the maximiser, not speed
# up obtention of initial conditions themselves; we'll still get RE ests
# from glmmTMB.

parseInits <- function(X, params, inds, inits.long){
  if(!inherits(X, "list"))
    stop(sQuote("inits"), " must be of class ", sQuote("list"), ".")
  
  poss.names <- c("D", "beta", "sigma", "gamma", "zeta")
  
  check.names <- sapply(names(X), `%in%`, poss.names)
  if(any(!check.names))
    stop(sQuote("inits"), " can only have elements: ", 
         paste0(sapply(poss.names, sQuote), collapse=', '), '.')
  
  # We'll proceed by 'overwriting' the initial conditions
  # obtained via glmmTMB/TimeVarCox. Leave those who aren't specified
  # in the control argument untouched.
  params.new <- params
  if(!is.null(X$D)){
    vD <- vech(X$D)
    wD <- grepl("^D\\[", names(params))
    if(length(vD)!=sum(wD)) # Check dim(D)
      stop("inits$D is improperly dimensioned.")
    params.new[wD] <- vD
  }
  if(!is.null(X$beta)){
    lb <- length(unlist(inds$R$beta))
    if(length(X$beta)!=lb)    # Check dim(beta)
      stop("inits$beta is improperly dimensioned.")
    beta.start <- max(grep("^D\\[", names(params))) + 1
    params.new[beta.start:(beta.start + lb - 1)] <- X$beta
  }
  if(!is.null(X$sigma)){
    if(!inherits(X$sigma, "list"))
      stop("Provide inits$sigma as a list of length K")
    if(length(X$sigma)!=length(inds$R$beta))
      stop("X$sigma must be of length K (set the unused elements to zero).")
    
    il.s <- inits.long$sigma.init
    il.s <- unlist(il.s)
    lS <- length(il.s)
    if(length(unlist(X$sigma))!=lS)
      stop("x$sigma improperly dimensioned.")
    
    to.replace <- names(il.s[il.s!=0])
    new.sigma <- unlist(X$sigma)
    new.sigma <- new.sigma[new.sigma != 0]
    
    # Assume they're in the correct order!
    params.new[which(names(params) == to.replace)] <- new.sigma
  }
  if(!is.null(X$gamma)){
    lg <- length(X$gamma)
    if(lg != length(inds$R$beta))
      stop("inits$gamma must be a vector of length K.")
    
    params.new[grepl("^gamma\\_", names(params))] <- X$gamma
  }
  if(!is.null(X$zeta)){
    lz <- length(X$zeta)
    wz <- which(grepl("^zeta\\_", names(params)))
    if(lz != length(wz))
      stop("inits$zeta improperly dimensioned vetor.")
    
    params.new[wz] <- X$zeta
  }
  
  params.new
}
