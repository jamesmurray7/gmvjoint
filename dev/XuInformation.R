ScoreOmega <- function(Perturbed.Omega, dmats, surv, sv, b, 
                       l0i, l0u, w, v, n, family, K, q, beta.inds, b.inds){
  D <- Perturbed.Omega$D
  beta <- c(Perturbed.Omega$beta)
  sigma <- Perturbed.Omega$sigma
  gamma <- c(Perturbed.Omega$gamma); gamma.rep <- rep(gamma, sapply(b.inds, length))
  zeta <- c(Perturbed.Omega$zeta)
  
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
  
  # Calculate b.hat, Sigma.hat at MLEs
  b.update <- Map(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    optim(b, joint_density, joint_density_ddb,
          Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
          beta_inds = beta.inds2, b_inds = b.inds2, K = K,
          method = 'BFGS', hessian = TRUE)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
  Fu = Fu, l0u = l0u)
  Sigma <- lapply(b.update, function(x) solve(x$hessian))
  b.hat <- lapply(b.update, function(x) x$par)
  SigmaSplit <- lapply(Sigma, function(x) lapply(b.inds, function(y) as.matrix(x[y,y])))
  bsplit <- lapply(b.hat, function(x) lapply(b.inds, function(y) x[y]))
  
  # Profile estimate for (gamma, zeta) at MLEs and bhat
  l0.hat <- lambdaUpdate_noprecalc(b.hat, Fu, SS, Sigma, gamma.rep, zeta, sv$nev, w, v)
  l0u.hat <- lapply(l0u, function(ll){
    l0.hat[1:length(ll)]
  })
  
  # Scores ------------------------------------------------------------------
  # The RE covariance matrix, D
  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = T), arr.ind = T)
  dimnames(vech.indices) <- NULL
  delta.D <- lapply(1:nrow(vech.indices), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1 # dD/dvech(d)_i i=1,...,length(vech(D))
    out
  })
  
  sDi <- function(i) {
    mapply(function(b, S) {
      EbbT <- S + tcrossprod(b)
      invDdD <- Dinv %*% delta.D[[i]]
      invDdDinvD <- invDdD %*% Dinv
      0.5 * sum(diag(invDdDinvD %*% EbbT)) - 0.5 * sum(diag(invDdD))
    },
    b = b.hat, S = Sigma,
    SIMPLIFY = TRUE)
  }
  
  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  
  tau <-  mapply(maketau, Z = Z, S = SigmaSplit, SIMPLIFY = F)
  
  Sb <- mapply(function(X, Y, Z, b, tau){
    c(Sbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K,
            FALSE, tau, w, v))
  }, X = X, Y = Y, Z = Z, b = bsplit, tau = tau, SIMPLIFY = F)
  
  # The dispersion parameter, \sigma
  Ss <- vector('list', K)
  funlist <- unlist(family)
  disps <- which(funlist %in% c('gaussian', 'genpois', 'Gamma'))
  for(j in disps){
    tauj <- lapply(tau, el, j)
    
    # Create score accordingly.
    if(funlist[j] == 'gaussian'){
      temp <- numeric(n)
      for(i in 1:n){
        rhs <- 0
        for(l in 1:length(w)){
          rhs <- rhs + w[l] * crossprod(Y[[i]][[j]] - X[[i]][[j]] %*% beta[beta.inds[[j]]] - Z[[i]][[j]] %*% bsplit[[i]][[j]] - v[l] * tauj[[i]])
        }
        temp[i] <- -m[[i]][j]/(2 * unlist(sigma)[j]) + 1/(2 * unlist(sigma)[j]^2) * rhs
      }
      Ss[[j]] <- temp
    }else if(funlist[j] == 'genpois'){
      Ss[[j]] <- mapply(function(b, X, Y, Z, tauj){
        phi_update(b[[j]], X[[j]], Y[[j]], Z[[j]], beta[beta.inds[[j]]], sigma[[j]], w, v, tauj)$Score
      }, b = bsplit, X = X, Y = Y, Z = Z, tauj = tauj, SIMPLIFY = F)
    }else if(funlist[j] == 'Gamma'){
      Ss[[j]] <- mapply(function(b, X, Y, Z, tauj){
        pracma::grad(E_shape.b, sigma[[j]], X = X[[j]], Y = Y[[j]], Z = Z[[j]], tauj = tauj, 
                     beta = beta[beta.inds[[j]]], b = b[[j]], w = w, v = v)
      }, b = bsplit, X = X, Y = Y, Z = Z, tauj = tauj, SIMPLIFY = F)
    }else{
      Ss[[j]] <- as.list(rep(NULL, n))
    }
  }
  
  # Survival parameters (\gamma, \zeta)
  Sgz <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta_cd(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = Sigma, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u.hat, Delta = Delta)
  
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

paramsToList <- function(params, q, P, K){
  Omega <- list()
  # D
  vD <- params[1:(q*(q+1)/2)]; upto <- (q*(q+1)/2) 
  Omega$D <- vech2mat(vD, q)
  # beta
  beta <- params[(upto + 1):(upto + P)]; upto <- upto + P
  Omega$beta <- c(beta)
  # sigma
  sigma <- c(params[(upto + 1):(upto + 1)],0); upto <- upto + 1
  Omega$sigma <- as.list(sigma)
  # gamma,zeta
  gamma <- params[(upto + 1):(upto + K)]; upto <- upto + K
  zeta <- params[(upto + 1):length(params)]
  Omega$gamma <- gamma; Omega$zeta <- zeta
  Omega
}

.ScoreWrapper <- function(params, dmats, surv, sv, b, 
                          l0i, l0u, w, v, n, family, K, q, beta.inds, b.inds){
  Omega.p <- paramsToList(params, q, sum(dmats$P), K)
  ScoreOmega(Omega.p, dmats, surv, sv, b, 
             l0i, l0u, w, v, n, family, K, q, beta.inds, b.inds)
}

Xu.Information <- function(params, dmats, surv, sv, b, 
                           l0i, l0u, w, v, n, family, K, q, beta.inds, b.inds,
                           method = "central", heps = 1e-4){
  numDiff(params, .ScoreWrapper, 
          dmats = dmats, surv = surv, sv = sv, b = b, l0i = l0i,
          l0u = l0u, w = w, v = v, n = n, family = family,
          K = K, q = q, beta.inds = beta.inds, b.inds = b.inds,
          method = method, heps = heps)
}

Itilde <- Xu.Information(params, dmats, surv, sv, b, 
               l0i, l0u, w, v, n, family, K, q, beta.inds, b.inds,
               method = 'central', heps = 1e-3)

Xu.SE <- setNames(sqrt(diag(solve(-Itilde))),
                  names(params))
