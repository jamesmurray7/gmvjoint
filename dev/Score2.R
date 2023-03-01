ll <- function(params, dmats, surv, sv, family){
  
  # Unpack the perturbed parameter vector
  D <- params[grepl('^D\\[', names(params))]
  lD <- length(D)
  D.tilde <- vech2mat(D, sv$q)
  # betas
  beta.tilde <- params[(lD + 1):(lD + sum(dmats$P))]
  # dispersions
  funlist <- unlist(family)
  sigmas <- params[(lD + sum(dmats$P) + 1):
                     (lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois')))]
  temp <- rep(NA, length(funlist))
  for(j in seq_along(funlist)){
    if(funlist[j] %in% c("gaussian", "Gamma", 'genpois')) temp[j] <- sigmas[j] else temp[j] <- 0
  }
  sigma.tilde <- as.list(temp)
  # gamma,zeta
  gz <- params[(1 + lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois'))):
                 length(params)] 
  gamma.tilde <- gz[1:length(dmats$np)]
  zeta.tilde <- gz[(1 + length(dmats$np)):length(gz)]
  gamma.tilde.rep <- rep(gamma.tilde, dmats$q)
  
  Omega.tilde <- list(D = D.tilde, beta = beta.tilde, sigma = sigma.tilde,
                      gamma = gamma.tilde, zeta = zeta.tilde)
  
  # Calculate {b.hat, lambda.hat} at the perturbed values
  # b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
  #     optim(b, joint_density, joint_density_ddb,
  #           Y = Y, X = X, Z = Z, beta = beta.tilde, D = D.tilde, sigma = sigma.tilde, family = family,
  #           Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.tilde.rep, zeta = zeta.tilde,
  #           beta_inds = beta.inds2, b_inds = b.inds2, K = K,
  #           method = 'L-BFGS-B', hessian = TRUE, lower = -Inf, upper = Inf)
  #   }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS,
  # Fu = sv$Fu, l0u = sv$l0u, SIMPLIFY = F)
  # b.hat.tilde <- lapply(b.update, el, 1)
  # Sigma.hat.tilde <- lapply(b.update, function(x) solve(x$hess))
    
  lambda.tilde <- sv$nev/rowSums(lambdaUpdate(sv$surv.times, sv$ft.mat, gamma.tilde, zeta.tilde, sv$S, 
                                              Sigma, b, w, v, b.inds2))
  l0u.tilde <- lapply(sv$l0u, function(ll){
    lambda.tilde[1:length(ll)]
  })
  l0i.tilde <- lambda.tilde[match(sv$Ti, sv$ft)] 
  l0i.tilde[is.na(l0i.tilde)] <- 0
  l0i.tilde <- as.list(l0i.tilde)
  
  ll.tilde <- mapply(function(Y, X, Z, b, S, SS, Fu, Fi, Delta, l0i, l0u){
    joint_density(b = b, Y = Y, X = X, Z = Z, beta = beta.tilde, D = D.tilde, sigma = sigma.tilde, family = family,
                  Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.tilde.rep, zeta = zeta.tilde,
                  beta_inds = beta.inds2, b_inds = b.inds2, K = K) * -1
  }, Y = dmats$Y, X = dmats$X, Z = dmats$Z, b = b, S = sv$S, SS = sv$SS, Fu = sv$Fu, Fi = sv$Fi, 
     Delta = surv$Delta, l0i = l0i.tilde, l0u = l0u.tilde)
  sum(ll.tilde)
}

# Takes absolutely ages (with b maximisation), but they look about right
Ioc <- pracma::hessian(ll, params, dmats = dmats, surv = surv, sv = sv, family = family, h = max(abs(params), 1) * 1e-7)
sqrt(diag(solve(-Ioc)))


.getScore <- function(Omega, dmats, surv, sv,
                                  Sigma, SigmaSplit, b, bsplit, 
                                  l0u, w, v, n, family, K, q, beta.inds, b.inds, beta.quad){
  # Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  sigma <- Omega$sigma
  gamma <- c(Omega$gamma)
  gamma.rep <- rep(gamma, sapply(b.inds, length))
  zeta <- c(Omega$zeta)
  
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
  
  #TEMP --> Consternation wrt if b should be re-maximised at Omega.tilde.
  # b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
  #     optim(b, joint_density, joint_density_ddb,
  #           Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
  #           Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
  #           beta_inds = beta.inds2, b_inds = b.inds2, K = K,
  #           method = 'BFGS', hessian = TRUE #lower = -Inf, upper = Inf,
  #           #, control = list(reltol=1e-4)
  #           )
  #   }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, Delta = surv$Delta, S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS,
  # Fu = sv$Fu, l0u = sv$l0u, SIMPLIFY = F)
  # b <- lapply(b.update, el, 1)
  # Sigma <- lapply(b.update, function(x) solve(x$hess))
  
  lambda.tilde <- sv$nev/rowSums(lambdaUpdate(sv$surv.times, sv$ft.mat, gamma, zeta, sv$S, 
                                              Sigma, b, w, v, b.inds2))
  l0u <- lapply(sv$l0u, function(ll){
    lambda.tilde[1:length(ll)]
  })
  
  # Scores ------------------------------------------------------------------
  # The RE covariance matrix, D
  Dinv <- solve(D)
  postmult <- diag(1, nrow = nrow(Dinv), ncol = ncol(Dinv))
  postmult[postmult == 0] <- 2 # off-diagonals have _twice_ the contribution!
  sD <- mapply(function(b, S){
    vech(t(0.5 * (Dinv %*% (S + tcrossprod(b)) %*% Dinv) - 0.5 * Dinv)) * vech(postmult)
  }, b = b, S = Sigma, SIMPLIFY = F)

  Sb <- mapply(function(X, Y, Z, b, tau){
    c(Sbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K,
            beta.quad, tau, w, v))
  }, X = X, Y = Y, Z = Z, b = bsplit, tau = list(0), SIMPLIFY = F)
  
  # The dispersion parameter, \sigma
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
          rhs <- rhs + w[l] * crossprod(Y[[i]][[j]] - X[[i]][[j]] %*% beta[beta.inds[[j]]] - Z[[i]][[j]] %*% bsplit[[i]][[j]] - v[l] * tau[[i]])
        }
        temp[i] <- -m[[i]][j]/(2 * unlist(sigma)[j]) + 1/(2 * unlist(sigma)[j]^2) * rhs
      }
      Ss[[j]] <- temp
    }else if(funlist[j] == 'genpois'){
      Ss[[j]] <- mapply(function(b, X, Y, Z, tau){
        phi_update(b[[j]], X[[j]], Y[[j]], Z[[j]], beta[beta.inds[[j]]], sigma[[j]], w, v, tau)$Score
      }, b = bsplit, X = X, Y = Y, Z = Z, tau = tau, SIMPLIFY = F)
    }else if(funlist[j] == 'Gamma'){
      Ss[[j]] <- mapply(function(b, X, Y, Z, tau){
        pracma::grad(E_shape.b, sigma[[j]], X = X[[j]], Y = Y[[j]], Z = Z[[j]], tau = tau, 
                     beta = beta[beta.inds[[j]]], b = b[[j]], w = w, v = v)
      }, b = bsplit, X = X, Y = Y, Z = Z, tau = tau, SIMPLIFY = F)
    }else{
      Ss[[j]] <- as.list(rep(NULL, n))
    }
  }
  
  # Survival parameters (\gamma, \zeta)
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, .Machine$double.eps^(1/3))
  }, b = b, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # Collate and form information --------------------------------------------
  Ss2 <- lapply(1:n, function(i){
    do.call(c, lapply(disps, function(j){
      Ss[[j]][[i]]
    }))
  })
  
  S <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, Sb, Ss, c(Sgz))
  }, sD = sD, Sb = Sb, Ss = Ss2, Sgz = Sgz)
  
  rowSums(S)
}


.Score <- function(params){
  
  # Unpack the perturbed parameter vector
  D <- params[grepl('^D\\[', names(params))]
  lD <- length(D)
  D.tilde <- vech2mat(D, sv$q)
  # betas
  beta.tilde <- params[(lD + 1):(lD + sum(dmats$P))]
  # dispersions
  funlist <- unlist(family)
  sigmas <- params[(lD + sum(dmats$P) + 1):
                     (lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois')))]
  temp <- rep(NA, length(funlist))
  for(j in seq_along(funlist)){
    if(funlist[j] %in% c("gaussian", "Gamma", 'genpois')) temp[j] <- sigmas[j] else temp[j] <- 0
  }
  sigma.tilde <- as.list(temp)
  # gamma,zeta
  gz <- params[(1 + lD + sum(dmats$P) + sum(unlist(family) %in% c("gaussian", "Gamma", 'genpois'))):
                 length(params)] 
  gamma.tilde <- gz[1:length(dmats$np)]
  zeta.tilde <- gz[(1 + length(dmats$np)):length(gz)]
  gamma.tilde.rep <- rep(gamma.tilde, dmats$q)
  
  Omega.tilde <- list(D = D.tilde, beta = beta.tilde, sigma = sigma.tilde,
                      gamma = gamma.tilde, zeta = zeta.tilde)
  return(.getScore(Omega.tilde, dmats, surv, sv, Sigma, SigmaSplit, b, bsplit, l0u, w, v
         , n, family, K, q, beta.inds, b.inds, beta.quad))

  
}


# "PFDS" (and appxs of it) from Xu 2014 -----------------------------------
qq <- function(x) setNames(sqrt(diag(solve(-x))), names(params))
# Forward differencing (fastest. pracma::Jacobian recommend heps = 1e-4).
J <- pracma::jacobian(.Score, params, heps = max(abs(params), 1) * 1e-4)
J <- 0.5 * (J + t(J))
qq(J)
J2c <- numDiff(params, .Score, method = 'central', heps = 1e-4)
J2f <- numDiff(params, .Score, method = 'forward', heps = 1e-3)
J2R <- cendiff2(params, .Score, method = 'Richardson', eps = 1e-4)

qq(J2c)
qq(J2f)
qq(J2R)
