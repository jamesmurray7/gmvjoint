# Observed data score ----
oScore <- function(params, dmats, surv, sv, family){
  # Unpack the parameter vector
  # D
  D <- matrix(0, nrow = sum(dmats$q), ncol = sum(dmats$q))
  D[lower.tri(D, T)] <- params[grepl('^D\\[', names(params))]
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  lD <- length(vech(D))
  
  # betas
  betas <- params[(lD + 1):(lD + sum(dmats$P))]
  # dispersions
  sigmas <- params[(lD + sum(dmats$P) + 1):
                     (lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois')))]
  sigmas <- as.list(sigmas)
  # gamma,zeta
  gz <- params[(1 + lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois'))):
                 length(params)] 
  gamma <- gz[1:length(dmats$np)]
  zeta <- gz[(1 + length(dmats$np)):length(gz)]
  gamma.rep <- rep(gamma, dmats$q)
  
  # b at perturbed Omega
  # b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
  #   optim(b, joint_density, joint_density_ddb,
  #         Y = Y, X = X, Z = Z, beta = betas, D = D, sigma = sigmas, family = family, 
  #         Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
  #         beta_inds = beta.inds2, b_inds = b.inds2, K = K,
  #         method = 'L-BFGS-B', hessian = TRUE, lower = -Inf, upper = Inf)
  # }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = l0i, SS = sv$SS,
  # Fu = sv$Fu, l0u = l0u, SIMPLIFY = F)
  # b <- lapply(b.update, el, 1)
  # Sigma <- lapply(b.update, function(x) solve(x$hess))
  
  # D ---------------------------------------------------------------------
  Dinv <- solve(D)
  postmult <- diag(1, nrow = nrow(Dinv), ncol = ncol(Dinv))
  postmult[postmult == 0] <- 2
  sD <- mapply(function(b, S){ # This one is correct, corresponds to what was done above, but quicker.
    vech(t(-(0.5 * Dinv - 0.5 * Dinv %*% (S + tcrossprod(b)) %*% Dinv))) * vech(postmult) # OBSERVED Score
  }, b = b, S = Sigma, SIMPLIFY = F)

  # beta --------------------------------------------------------------------
  fnb <- function(bt, Y, X, Z, b, S, w, v){
    tau <- sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]])))
    eta <- X[[1]] %*% bt + Z[[1]] %*% b
    out <- matrix(0, nr = length(Y[[1]]), nc = length(w))
    for(l in 1:length(w)){
      eta.l <- eta + tau * v[l]
      out[,l] <- w[l] * dnorm(Y[[1]], eta.l, sigmas[[1]], log = T)
    }
    sum(rowSums(out))
  }
  
  Sb <- mapply(function(Y, X, Z, b, S){
    pracma::grad(fnb, betas, Y = Y, X = X, Z = Z, b = b, S = S, w = w, v = v)
  }, Y = Y, X = X, Z = Z, b = b, S = Sigma, SIMPLIFY = F)
  
  
  # sigma -------------------------------------------------------------------
  Ss <- vector('list', K)
  funlist <- unlist(family)
  disps <- which(funlist %in% c('gaussian', 'genpois', 'Gamma'))
  for(j in disps){
    if(beta.quad)
      tau <- lapply(tau, el, j)
    else
      tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z[[j]] %*% S[[j]], Z[[j]])))), Z = Z, S = SigmaSplit)
    
    # Create score accordingly.
    if(funlist[j] == 'gaussian'){
      temp <- numeric(n)
      for(i in 1:n){
        rhs <- 0
        for(l in 1:length(w)){
          rhs <- rhs + w[l] * crossprod(Y[[i]][[j]] - X[[i]][[j]] %*% betas[beta.inds[[j]]] - Z[[i]][[j]] %*% bsplit[[i]][[j]] - v[l] * tau[[i]])
        }
        temp[i] <- -m[[i]][j]/(2 * unlist(sigmas)[j]) + 1/(2 * unlist(sigmas)[j]^2) * rhs
      }
      Ss[[j]] <- temp
    }else if(funlist[j] == 'genpois'){
      Ss[[j]] <- mapply(function(b, X, Y, Z, tau){
        phi_update(b[[j]], X[[j]], Y[[j]], Z[[j]], betas[beta.inds[[j]]], sigmas[[j]], w, v, tau)$Score
      }, b = bsplit, X = X, Y = Y, Z = Z, tau = tau, SIMPLIFY = F)
    }else if(funlist[j] == 'Gamma'){
      Ss[[j]] <- mapply(function(b, X, Y, Z, tau){
        pracma::grad(E_shape.b, sigmas[[j]], X = X[[j]], Y = Y[[j]], Z = Z[[j]], tau = tau, 
                     beta = betas[beta.inds[[j]]], b = b[[j]], w = w, v = v)
      }, b = bsplit, X = X, Y = Y, Z = Z, tau = tau, SIMPLIFY = F)
    }else{
      Ss[[j]] <- as.list(rep(NULL, n))
    }
  }
  
  # (gamma, zeta) -----------------------------------------------------------
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Ss2 <- lapply(1:n, function(i){
    do.call(c, lapply(disps, function(j){
      Ss[[j]][[i]]
    }))
  })
  
  Scores <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, Sb, Ss, c(Sgz))
  }, sD = sD, Sb = Sb, Ss = Ss2, Sgz = Sgz)
  
  Scores
}

rSoS <- function(Omega) rowSums(oScore(params = Omega, dmats, surv, sv, family))
rSoS(params)
temp <- cendiff(params, rSoS)
sqrt(diag(solve(-temp))) # This is 4.7.1 Bakers' method (? But M&K cite Meilijson (1989)).

Score3 <- function(params, dmats, surv, sv, Sigma, SigmaSplit, b, bsplit, l0u, w, v, n, family, 
                   K, q, beta.inds, b.inds){
  D <- matrix(0, nrow = sum(dmats$q), ncol = sum(dmats$q))
  D[lower.tri(D, T)] <- params[grepl('^D\\[', names(params))]
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  lD <- length(vech(D))
  
  # betas
  betas <- params[(lD + 1):(lD + sum(dmats$P))]
  # dispersions
  sigmas <- params[(lD + sum(dmats$P) + 1):
                     (lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois')))]
  sigmas <- as.list(sigmas)
  # gamma,zeta
  gz <- params[(1 + lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois'))):
                 length(params)] 
  gamma <- gz[1:length(dmats$np)]
  zeta <- gz[(1 + length(dmats$np)):length(gz)]
  
  # Extract data objects ----
  # Longitudinal //
  Z <- dmats$Z
  X <- dmats$X
  Y <- dmats$Y
  m <- lapply(Y, function(y) sapply(y, length))
  
  beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) # Indexed for C++ use.
  
  # Survival //
  S <- sv$S
  SS <- sv$SS
  Fi <- sv$Fi
  Fu <- sv$Fu
  Delta <- surv$Delta
  

  # D -----------------------------------------------------------------------
  Dinv <- solve(D)
  postmult <- diag(1, nrow = nrow(Dinv), ncol = ncol(Dinv))
  postmult[postmult == 0] <- 2 # off-diagonals have twice the contribution!
  # Observed
  sDobs <- mapply(function(b, S){ 
    vech(t(-(0.5 * Dinv - 0.5 * Dinv %*% (S + tcrossprod(b)) %*% Dinv))) * vech(postmult) # OBSERVED Score
  }, b = b, S = Sigma, SIMPLIFY = F)
  # Complete
  sDcom <- mapply(function(b){ 
    vech(t(-(0.5 * Dinv - 0.5 * Dinv %*% (tcrossprod(b)) %*% Dinv))) * vech(postmult) # COMPLETE Score
  }, b = b, SIMPLIFY = F)
  
  # Checking with pracma::grad
  temp.D <- function(vD, b, S, obs = T){
    D <- diag(2)
    D[lower.tri(D, T)] <- vD
    D[upper.tri(D)] <- t(D)[upper.tri(D)]
    
    if(obs)
      return(-0.5 * log(det(D)) - 0.5 * sum(diag(solve(D) %*% (S + tcrossprod(b)))))
    else
      return(-0.5 * log(det(D)) - 0.5 * crossprod(b, solve(D) %*% b))
  }
  # pracma::grad(temp.D, vD, b = b[[1]], S = Sigma[[1]], obs = T) # Check with sD[[1]]
  
  # beta --------------------------------------------------------------------
  fnb <- function(bt, Y, X, Z, b, S, w, v){
    tau <- sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]])))
    eta <- X[[1]] %*% bt + Z[[1]] %*% b
    out <- matrix(0, nr = length(Y[[1]]), nc = length(w))
    for(l in 1:length(w)){
      eta.l <- eta + tau * v[l]
      out[,l] <- w[l] * dnorm(Y[[1]], eta.l, sigmas[[1]], log = T)
    }
    sum(rowSums(out))
  }
  Sbobs <- mapply(function(Y, X, Z, b, S){
    pracma::grad(fnb, betas, Y = Y, X = X, Z = Z, b = b, S = S, w = w, v = v)
  }, Y = Y, X = X, Z = Z, b = b, S = Sigma, SIMPLIFY = F)
  Sbcomp <- mapply(function(Y, X, Z, b, S){
    pracma::grad(fnb, betas, Y = Y, X = X, Z = Z, b = b, S = S, w = c(1), v = c(0))
  }, Y = Y, X = X, Z = Z, b = b, S = Sigma, SIMPLIFY = F)
  
  # The dispersion parameter, \sigma
  Ssobs <- Sscomp <- vector('list', n)
  tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z[[1]] %*% S, Z[[1]])))), Z = Z, S = Sigma)
  for(i in 1:n){
    comprhs <- crossprod(Y[[i]][[1]] - X[[i]][[1]] %*% betas - Z[[i]][[1]] %*% b[[i]])
    Sscomp[[i]] <- -m[[i]][1]/(2 * unlist(sigmas)[1]) + 1/(2 * unlist(sigmas)[1]^2) * comprhs
    obsrhs <- 0
    for(l in 1:length(w)){
      obsrhs <- obsrhs + w[l] * crossprod(Y[[i]][[1]] - X[[i]][[1]] %*% betas - Z[[i]][[1]] %*% b[[i]] - tau[[i]] * v[l])
    }
    Ssobs[[i]] <- -m[[i]][1]/(2 * unlist(sigmas)[1]) + 1/(2 * unlist(sigmas)[1]^2) * obsrhs
  }
  

  # gamma, zeta -------------------------------------------------------------
  Sgzobs <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Sgzcomp <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w = c(1), v = c(0), b.inds2, K, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Imissing <- mapply(function(sdc, sdo, sbc, sbo, ssc, sso, sgzc, sgzo){
    com <- c(sdc, sbc, ssc, sgzc)
    obs <- c(sdo, sbo, sso, sgzo)
    tcrossprod(obs) - tcrossprod(com)
  }, sdc = sDcom, sdo = sDobs, sbc = Sbcomp, sbo = Sbobs,
     ssc = Sscomp, sso = Ssobs,
     sgzc = Sgzcomp, sgzo = Sgzobs, SIMPLIFY = F)
  Im <- Reduce('+', Imissing)
  
  
  # Complete information ----------------------------------------------------
  .jd <- function(params, b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz){
    D <- matrix(0, nrow = sum(dmats$q), ncol = sum(dmats$q))
    D[lower.tri(D, T)] <- params[grepl('^D\\[', names(params))]
    D[upper.tri(D)] <- t(D)[upper.tri(D)]
    lD <- length(vech(D))
    
    # betas
    betas <- params[(lD + 1):(lD + sum(dmats$P))]
    # dispersions
    sigmas <- params[(lD + sum(dmats$P) + 1):
                       (lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois')))]
    sigmas <- as.list(sigmas)
    # gamma,zeta
    gz <- params[(1 + lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois'))):
                   length(params)] 
    gamma <- gz[1:length(dmats$np)]
    zeta <- gz[(1 + length(dmats$np)):length(gz)]
    gamma.rep <- rep(gamma, dmats$q)
    joint_density(b, Y, X, Z, betas, D, sigmas, family, Delta, S, Fi, l0i, SS,
                  Fu, haz, gamma.rep, zeta, beta.inds2, b.inds2, K) * -1
  }
  Ic <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz){
    -1 * pracma::hessian(.jd, params,
                    b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, 
                    l0i = l0i, SS = SS, Fu = Fu, haz = haz)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, 
  l0i = l0i, SS = SS, Fu = Fu, haz = l0u, SIMPLIFY = F)
  
  sqrt(diag(solve(Reduce('+', Ic) - Im)))
  
  # Louis'
  Sc <- colSums(do.call(rbind,mapply(function(sdc, sbc, ssc, sgzc){
    c(sdc, sbc, ssc, sgzc)
  }, sdc = sDobs, sbc = Sbobs, ssc = Ssobs, sgzc = Sgzobs, SIMPLIFY = F)))
  
  sqrt(diag(solve(Reduce('+', Ic) - tcrossprod(Sc))))
  
}


# Louis -------------------------------------------------------------------
lse <- function(params, b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz){
  D <- matrix(0, nrow = sum(dmats$q), ncol = sum(dmats$q))
  D[lower.tri(D, T)] <- params[grepl('^D\\[', names(params))]
  D[upper.tri(D)] <- t(D)[upper.tri(D)]
  lD <- length(vech(D))
  
  # betas
  betas <- params[(lD + 1):(lD + sum(dmats$P))]
  # dispersions
  sigmas <- params[(lD + sum(dmats$P) + 1):
                     (lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois')))]
  sigmas <- as.list(sigmas)
  # gamma,zeta
  gz <- params[(1 + lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois'))):
                 length(params)] 
  gamma <- gz[1:length(dmats$np)]
  zeta <- gz[(1 + length(dmats$np)):length(gz)]
  gamma.rep <- rep(gamma, dmats$q)
  joint_density(b, Y, X, Z, betas, D, sigmas, family, Delta, S, Fi, l0i, SS,
                Fu, haz, gamma.rep, zeta, beta.inds2, b.inds2, K) * -1
}

SH <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz){
  out <- vector('list', 2)
  out[[1]] <- pracma::grad(lse, params,
                           b = b, Y = Y, X = X, Z = Z, Delta = Delta,
                           S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = haz)
  out[[2]] <- pracma::hessian(lse, params,
                           b = b, Y = Y, X = X, Z = Z, Delta = Delta,
                           S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = haz)
  out
}, b = b, Y = Y, X = X, Z = Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i,
   SS = sv$SS, Fu = sv$Fu, haz = sv$l0u, SIMPLIFY = F)

sqrt(diag(solve(-1 * (Reduce('+', lapply(SH, el , 2)) - tcrossprod(colSums(do.call(rbind, lapply(SH, el, 1))))))))
