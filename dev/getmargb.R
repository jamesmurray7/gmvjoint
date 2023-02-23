# Obtain f(b|Y,T,Delta;Omega) at given joint model fit by Metropolis scheme.
get.marg.b <- function(fit, burnin = 500L, N = 3500L, tune = 2.){
  if(!inherits(fit, "joint")) stop("Only usable with object of class 'joint'.")
  if(is.null(fit$dmats)) stop("Need dmats.")
  # Unpack dmats
  M <- fit$ModelInfo
  dm <- fit$dmats$long
  sv <- fit$dmats$surv
  surv <- fit$dmats$ph
  q <- sv$q; K <- length(M$family)
  
  # Unpack parameter estimates
  D <- fit$coeffs$D
  beta <- fit$coeffs$beta
  sigma <- fit$coeffs$sigma
  gamma <- fit$coeffs$gamma
  gamma.rep <- rep(gamma, sapply(M$inds$b, length))
  zeta <- fit$coeffs$zeta
  # Model matrices
  l0i <- sv$l0i; l0u <- sv$l0u
  Del <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi
  S <- sv$S; SS <- sv$SS
  X <- dm$X; Y <- dm$Y; Z <- dm$Z
  # Other
  beta.inds <- lapply(M$inds$beta, function(x) x-1)
  b.inds <- lapply(M$inds$b, function(x) x-1)
  ff <- M$family
  nobs <- sapply(lapply(Y, el), length)
  
  # Inits + metropolis scheme inputs
  b <- lapply(1:M$n, function(i) c(fit$REs[i,,drop=F]))
  Sigmas <- lapply(1:M$n, function(i) attr(fit$REs, 'vcov')[i,])
  Sigmas <- lapply(Sigmas, vech2mat, q)
  iters <- burnin + N
  out <- vector('list', M$n)
  accepts <- numeric(M$n)
  pb <- utils::txtProgressBar(max = M$n, style = 3)
  start.time <- proc.time()[3]
  for(i in 1:M$n){
    store <- matrix(0, nrow = iters, ncol = q)
    b.current <- b[[i]]
    accept <- numeric(iters)
    # Metropolis scheme
    for(j in 1:iters){
      b.prop <- MASS::mvrnorm(n = 1, mu = b.current, Sigma = Sigmas[[i]] * tune)
      fyTb.current <- exp(-joint_density(b.current, Y = Y[[i]], X = X[[i]], Z = Z[[i]],
                                         beta = beta, D = D, sigma = sigma, family = ff, 
                                         Delta = Del[[i]], S = S[[i]], Fi = Fi[[i]], l0i = l0i[[i]],
                                         SS = SS[[i]], Fu = Fu[[i]],  haz = l0u[[i]], gamma_rep = gamma.rep,
                                         zeta = zeta, beta_inds = beta.inds, b_inds = b.inds, K = K))
      fyTb.proposl <- exp(-joint_density(b.prop, Y = Y[[i]], X = X[[i]], Z = Z[[i]],
                                         beta = beta, D = D, sigma = sigma, family = ff, 
                                         Delta = Del[[i]], S = S[[i]], Fi = Fi[[i]], l0i = l0i[[i]],
                                         SS = SS[[i]], Fu = Fu[[i]],  haz = l0u[[i]], gamma_rep = gamma.rep,
                                         zeta = zeta, beta_inds = beta.inds, b_inds = b.inds, K = K))
      P <- min(fyTb.proposl/fyTb.current, 1)
      U <- runif(1)
      if(U <= P){
        b.current <- b.prop
        accept[j] <- 1
      }
      store[j,] <- b.current
    }
    utils::setTxtProgressBar(pb, i)
    out[[i]] <- store[-(1:burnin),]
    accept <- accept[-(1:burnin)]
    accepts[i] <- sum(accept)/N
  }
  end.time <- proc.time()[3]
  close(pb)
  out <- list(walks = out, acceptance = accepts, M = M, bhats = do.call(rbind, b),
              q = q, K = K, qnames = colnames(fit$coeffs$D),
              burnin = burnin, N = N, tune = tune, nobs = nobs,
              elapsed.time = end.time-start.time)
  class(out) <- 'marginal.b.joint'
  out
}

#' @keywords internal
#' @method print marginal.b.joint
print.marginal.b.joint <- function(x, ...){
  if(!inherits(x, 'marginal.b.joint')) stop("x must be a 'marginal.b.joint' object.")
  M <- x$M
  # Data information
  cat("Marginal distribution f(b_i|Y_i,T_i,Delta_i;Omega)\n")
  cat("------------------------------\n")
  cat(sprintf("Number of subjects: %d\n", M$n))
  cat(sprintf("Number of events: %d (%.2f%%)\n", M$nev, 100 * M$nev/M$n))
  cat(sprintf("Number of responses: %d, dimension of random effects: %d\n", x$K, x$q))
  nobs <- x$nobs
  qnobs <- quantile(nobs)
  if(qnobs[3] == qnobs[2] && qnobs[3] == qnobs[4])
    cat(sprintf("Median profile length: %d\n", qnobs[3]))
  else
    cat(sprintf("Median [IQR] profile length: %d [%d, %d]\n", qnobs[3], qnobs[2], qnobs[4]))
  
  # MC Information
  cat("\n")
  accepts <- x$acceptance
  walks <- do.call(rbind, x$walks)
  colnames(walks) <- x$qnames
  cat("Metropolis Algorithm information:\n")
  cat(sprintf("Burn-in: %d, MC iterations post burn-in: %d\n", x$burnin, x$N))
  cat(sprintf("Median acceptance rate: %.2f%%\n", median(accepts) * 100))
  
  # Print summary of the things 
  cat("\n")
  cat("Summary: \n")
  summ <- apply(walks, 2, summary)
  summ <- rbind(summ[1:4,], Emp.SD = apply(walks, 2, sd), summ[5:6,])
  print(round(summ, 4))
  
  invisible(x)
}

#' @keywords internal
#' @method plot marginal.b.joint
plot.marginal.b.joint <- function(x, D = NULL, nrow = NULL, ncol = NULL, title = NULL, ...){
  if(!inherits(x, 'marginal.b.joint')) stop("x must be a 'marginal.b.joint' object.")
  .par <- par(no.readonly = T)
  walks <- do.call(rbind, x$walks)
  num.plots <- x$q # Determine number of plots
  if(is.null(nrow) & !is.null(ncol) | !is.null(nrow) & is.null(ncol)) 
    stop("Define either both 'nrow' and 'ncol', or leave both unspecified.")
  if(is.null(nrow) & is.null(ncol)){ # Plot as if intercept and slope by default.
    nrow <- max(cumsum(unlist(x$M$ind$b)%%2))
    ncol <- 2
  }
  par(mfrow = c(nrow, ncol))
  for(j in 1:num.plots){
    dj <- density(walks[,j])
    if(!is.null(D)){
      seqx <- seq(min(dj$x), max(dj$x), length.out = 1e3)
      dn <- dnorm(seqx, mean = 0, sd = sqrt(D[j,j]))
      ylims <- c(0, max(max(dj$y), max(dn)) + .1)
    }else{
      ylims <- c(min(dj$y), max(dj$y) + 0.1)
    }
    plot(dj, xlab = '', main = x$qnames[j],
         ylim = ylims)
    if(!is.null(D))
      lines(seqx, dn, lty = 3, col = 'red2')
  }
  if(!is.null(title))
    mtext(title, side - 3, line = -2, outer = TRUE)
  on.exit(par(.par))
}

get.marg.b.cpp <- function(fit, burnin = 500L, N = 3500L, tune = 2.){
  if(!inherits(fit, "joint")) stop("Only usable with object of class 'joint'.")
  if(is.null(fit$dmats)) stop("Need dmats.")
  # Unpack dmats
  M <- fit$ModelInfo
  dm <- fit$dmats$long
  sv <- fit$dmats$surv
  surv <- fit$dmats$ph
  q <- sv$q; K <- length(M$family)
  
  # Unpack parameter estimates
  Omega <- fit$coeffs
  gamma.rep <- rep(Omega$gamma, sapply(M$inds$b, length))
  
  # Model matrices
  l0i <- sv$l0i; l0u <- sv$l0u
  Del <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi
  S <- sv$S; SS <- sv$SS
  X <- dm$X; Y <- dm$Y; Z <- dm$Z
  # Other
  beta.inds <- lapply(M$inds$beta, function(x) x-1)
  b.inds <- lapply(M$inds$b, function(x) x-1)
  ff <- M$family
  nobs <- sapply(lapply(Y, el), length)
  
  # Inits + metropolis scheme inputs
  b <- lapply(1:M$n, function(i) c(fit$REs[i,,drop=F]))
  Sigmas <- lapply(1:M$n, function(i) attr(fit$REs, 'vcov')[i,])
  Sigmas <- lapply(Sigmas, vech2mat, q)
  iters <- burnin + N
  out <- vector('list', M$n)
  pb <- utils::txtProgressBar(max = M$n, style = 3)
  start.time <- proc.time()[3]
  for(i in 1:M$n){
    this <- metropolis(b[[i]], Omega, Y[[i]], X[[i]], Z[[i]], ff, 
               Del[[i]], S[[i]], Fi[[i]], l0i[[i]], SS[[i]], 
               Fu[[i]], l0u[[i]], gamma.rep, beta.inds, b.inds, 
               K, q, burnin, N, Sigmas[[i]], tune)
    utils::setTxtProgressBar(pb, i)
    out[[i]] <- t(this$walks)
    accepts[i] <- this$AcceptanceRate
  }
  end.time <- proc.time()[3]
  close(pb)
  out <- list(walks = out, acceptance = accepts, M = M, bhats = do.call(rbind, b),
              q = q, K = K, qnames = colnames(fit$coeffs$D),
              burnin = burnin, N = N, tune = tune, nobs = nobs,
              elapsed.time = end.time-start.time)
  class(out) <- 'marginal.b.joint'
  out
}
