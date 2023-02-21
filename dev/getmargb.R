# Obtain f(b|Y,T,Delta;Omega) at given joint model fit by Metropolis scheme.
get.marg.b <- function(fit, burnin, N, tune){
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
  
  # Inits + metropolis scheme inputs
  b <- lapply(1:M$n, function(i) c(fit$REs[i,,drop=F]))
  iters <- burnin + N
  out <- vector('list', M$n)
  accepts <- numeric(M$n)
  pb <- utils::txtProgressBar(max = M$n, style = 3)
  for(i in 1:M$n){
    store <- matrix(0, nrow = N, ncol = q)
    b.current <- b[[i]]
    # Burnin phase (if req.)
    if(burnin > 0){
      for(j in 1:burnin){
        b.prop <- MASS::mvrnorm(n = 1, mu = b.current, Sigma = D)
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
        P <- fyTb.proposl/fyTb.current
        U <- runif(1)
        if(U <= P){
          b.current <- b.prop
        }
      }
    }
    # Metropolis
    accept <- 0
    for(j in 1:N){
      b.prop <- MASS::mvrnorm(n = 1, mu = b.current, Sigma = D)
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
      P <- fyTb.proposl/fyTb.current
      U <- runif(1)
      if(U <= P){
        b.current <- b.prop
        accept <- accept + 1
      }
      store[j,] <- b.current
    }
    utils::setTxtProgressBar(pb, i)
    out[[i]] <- store
    accepts[i] <- accept/N
  }
  close(pb)
  list(walks = out, acceptance = accepts)
}

# Get `fit` from joint.Rd for example.
# aad <- get.marg.b(fit, 500, 5000, F)
# par(mfrow=c(2,2)); for(j in 1:4) plot(density(do.call(rbind, aad$walks)[,j]),main=bquote(b[.(j-1)]));par(mfrow=c(1,1))
