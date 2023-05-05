# Version of simData where the user can set baseline hazard to be 
# either Gompertz, Weibull or Exponential (with values of theta then not/being used)!
simData <- function(n = 250, ntms = 10, fup = 5, 
                    family = list('gaussian', 'gaussian'), 
                    sigma = c(0.16, 0.16),
                    beta = rbind(c(1, 0.10, 0.33, -0.50), c(1, 0.10, 0.33, -0.50)), D = NULL,
                    gamma = c(0.5, -0.5), zeta = c(0.05, -0.30),
                    basehaz = c("Gompertz","Weibull","Exp"), theta = c(-4, 0.2), cens.rate = exp(-3.5),
                    unif.times = TRUE,
                    dof = Inf,
                    random.formula = NULL,
                    return.ranefs = FALSE){
  
  # Checks --------------
  # Check family is valid option
  family <- lapply(family, function(f){
    if("function"%in%class(f)) f <- f()$family
    if(!f%in%c("gaussian", "binomial", "poisson", "genpois", "Gamma")) 
      stop('Family must be one of "gaussian", "binomial", "poisson", "genpois" or "Gamma".')
    f
  })
  # Check dispersion supplied if neg. binom.
  family <- lapply(family, function(f) family <- match.arg(f, c("gaussian", "binomial", "poisson", "genpois", "Gamma"), several.ok = F))
  
  # Checks wrt the fixed effects and dispersion parameters
  funlist <- unlist(family)
  num.gp <- length(which(funlist == 'genpois'))
  num.ga <- length(which(funlist == 'gaussian'))
  num.Ga <- length(which(funlist == 'Gamma'))
  
  # Fixed effects
  if(!"matrix"%in%class(beta) & length(family) == 1L) beta <- matrix(beta,nr=1) # Assume the provided numeric is the 1-dim case
  K <- nrow(beta)
  if(K != length(funlist)) stop('Incorrect number of families provided and/or beta terms.')
  if(K != length(gamma)) stop('Incorrect number of association parameters provided wrt beta terms.')
  if(is.null(D)) D <- diag(K * 2) else D <- D
  if(!is.null(random.formula)){
    if(length(random.formula) != K) stop('Provided random formulas must be a list of length K.')
  }else{
    if((K*2) != ncol(D)) stop('Incorrect dimension on supplied covariance matrix D, do you need to supply random.formula?')
  }
  
  # Check covariance matrix D is positive semi-definite.
  if(any(eigen(D)$value < 0) || det(D) <= 0) stop('D must be positive semi-definite')
  
  # Necessary parameters & data generation ----
  if(unif.times){
    time <- seq(0, fup, length.out = ntms)
    tau <- fup + 0.1
  }else{
    time <- c(replicate(n, c(0, sort(runif(ntms - 1, max = fup)))))
    tau <- fup
  }

  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = if(unif.times) rep(time, n) else time,
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices, used across ALL responses ----
  X <- model.matrix(~ time + cont + bin, data = df)
  if(!is.null(random.formula)){
    Z <- lapply(random.formula, function(x) model.matrix(x, data = df))
  }else{
    Z <- replicate(K, model.matrix(~ 1 + time, data = df), simplify = F) # assume all intslopes.
  }

  # Random effects specification
  if(is.infinite(dof)) 
    .simRE <- function(n, mu, Sigma) MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  else
    .simRE <- function(n, mu, Sigma) mvtnorm::rmvt(n = n, sigma = Sigma, delta = mu, df = dof)
  
  # Linear predictor & response generation ----
  if(!is.null(random.formula)){
    q <- ncol(do.call(cbind, lapply(Z, head)))
    # b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
    b <- .simRE(n, rep(0, q), D)
    b.inds <- split(1:q, do.call(c, sapply(1:K, function(k) rep(k, ncol(Z[[k]])), simplify = F)))
  }else{
    # b <- MASS::mvrnorm(n, mu = rep(0, K * 2), Sigma = D)
    b <- .simRE(n, rep(0, K * 2), D)
    b.inds <- split(1:(2*K), rep(1:K, each = 2)) 
  }
  
  Y <- sapply(1:K, function(k){
    f <- family[[k]]; Zk <- Z[[k]]
    betak <- beta[k,,drop=T]
    etak <- X %*% betak + rowSums(Zk * b[df$id, b.inds[[k]]])
    switch(f, 
           gaussian = Y <- rnorm(n * ntms, mean = etak, sd = sqrt(sigma[k])),
           binomial = Y <- rbinom(n * ntms, 1, plogis(etak)),
           poisson = Y <- rpois(n * ntms, exp(etak)),
           Gamma = Y <- rgamma(n * ntms, shape = sigma[k], scale = exp(etak)/sigma[k]),
           genpois = {
             maxtry <- 5
             try <- 1
             e <- T
             while(e){ # genpois can be quite sensitive in simulation part; repeat sim up to 5 times.
               Y <- tryCatch(rgenpois(exp(etak), sigma[k]),
                             error = function(e) NA)
               e <- any(is.na(Y))
               try <- try + 1
               if(try > maxtry) stop("Issues creating genpois response with beta values ", c(betak), 
                                     " and mean(b)", c(round(mean(b[df$id, b.inds[[k]]]), 3)))
             }
             Y
           })
  })
  colnames(Y) <- paste0('Y.', 1:K)
  
  df <- cbind(df, Y)
  
  # Survival ----
  # theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% zeta
  U <- runif(n)
  if(!is.null(random.formula)){
    zz <- lapply(Z, colnames)
    ints <- which(grepl('\\(Intercept\\)', do.call(c, zz))); slopes <- which(grepl('time', do.call(c, zz)))
    if(any(!grepl('\\(Intercept\\)|time', do.call(c, zz)))) message('Warning: Only intercept-only or intercept-ands-slope random.formulas allowed.')
    # Which of the K responses have intercept/slope
    ints.gamma <- do.call(c, lapply(1:K, function(k){
      if(any(grepl('\\(Intercept\\)', colnames(Z[[k]])))) return(k) else return(NULL)
    }))
    slopes.gamma <- do.call(c, lapply(1:K, function(k){
      if(any(grepl('time', colnames(Z[[k]])))) return(k) else return(NULL)
    }))
  }else{
    ints <- seq(1,2*K,by=2); slopes <- seq(2, 2*K, by = 2) # intslope on all.
    ints.gamma <- slopes.gamma <- 1:K    # Each K has both intercept and slope.
  }
  
  b0 <- b[, ints, drop = F]; b1 <- b[, slopes, drop = F]
  Q <- b1 %*% gamma[slopes.gamma]
  P <- Keta + b0 %*% gamma[ints.gamma]
  nu <- theta[1]; alpha <- theta[2] # Scale, shape
  # Generate survival times (Austin et al 2012)
  basehaz <- match.arg(basehaz)
  t <- switch(basehaz,
              Exp = {
                suppressWarnings(log(1-(Q * log(U))/exp(nu + P))/Q)
              },
              Weibull = {
                suppressWarnings(
                  (log(1-((1+alpha) * log(U))/(Q*alpha*exp(nu+P)))/Q)^(1/(1+alpha))
                )
              },
              Gompertz = {
                suppressWarnings(
                  log(1-((Q+alpha)*log(U))/exp(nu+P))/(Q+alpha)
                )
              })
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate cenors times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  # Output Dataframes ----
  surv.data <- data.frame(id = 1:n, survtime, status)
  long.data <- df
  
  out.data <- merge(df, surv.data, by = 'id')
  out.data <- out.data[out.data$time < out.data$survtime, ]
  cat(sprintf("%.2f%% failure rate.\n",round(100 * sum(surv.data$status)/n)))
  
  out <- list(data =  out.data, 
              surv.data =  out.data[!duplicated(out.data[,'id']), c('id', 'survtime', 'status', 'cont', 'bin')])
  if(return.ranefs) out$ranefs <- b
  out
}
