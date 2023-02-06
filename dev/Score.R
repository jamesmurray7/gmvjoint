Score <- function(params, dmats, surv, sv, family){ 
  
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
  

  # Score for D -------------------------------------------------------------
  
  # Scores ------------------------------------------------------------------
  # The RE covariance matrix, D
  Dinv <- solve(D)
  vech.indices <- which(lower.tri(D, diag = T), arr.ind = T)
  dimnames(vech.indices) <- NULL
  delta.D <- lapply(1:nrow(vech.indices), function(d){
    out <- matrix(0, nrow(D), ncol(D))
    ind <- vech.indices[d, 2:1]
    out[ind[1], ind[2]] <- out[ind[2], ind[1]] <- 1 # dD/dvech(d)_i
    out
  })
  
  lhs <- sapply(delta.D, function(d) {
    -0.5 * sum(diag(Dinv %*% d))
  })
  
  sDi <- function(i) {
    mapply(function(b, S) {
      out <- 0.5 * tcrossprod(b) %*% (Dinv %*% delta.D[[i]] %*% Dinv)   
      lhs[i] + sum(diag(out))
    },
    b = b, S = Sigma,
    SIMPLIFY = T)
  }
  
  sD <- colMeans(sapply(1:nrow(vech.indices), sDi))
  
  T0 <- Dinv
  T1 <- .5 * T0
  
  test1 <- mapply(function(b){
    t2 <- T0 %*% b
    a <- t(T1) + T1 - 0.5 * t(T0) %*% b %*% t(t2) - 0.5 * t2 %*% crossprod(b, T0)
    -0.5 * vech(a)
  },b = b, SIMPLIFY = F)
  
  test1b <- mapply(function(b){
    -(t(0.5 * T0) - 0.5 * t(T0) %*% b %*% t(T0 %*% b))
  }, b = b, SIMPLIFY = F)
  
  ScoreD <- function(vD, b){
    Dmat <- matrix(0, nrow = sum(dmats$q), ncol = sum(dmats$q))
    Dmat[lower.tri(Dmat, T)] <- vD
    Dmat[upper.tri(Dmat)] <- t(Dmat)[upper.tri(Dmat)]
    
    
    mvtnorm::dmvnorm(b, sigma = Dmat, log = T)
  }
  
  test2 <- lapply(1:n, function(i) pracma::grad(ScoreD, vD, b = b[[i]]))
  
  
  # The fixed effects, \beta 
  if(beta.quad){
    tau = mapply(maketau, Z = Z, S = SigmaSplit, SIMPLIFY = F)
  }else{
    tau = list(0)
  }
  
  Sb <- mapply(function(X, Y, Z, b, tau){
    c(Sbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K,
            beta.quad, tau, w, v))
  }, X = X, Y = Y, Z = Z, b = bsplit, tau = tau, SIMPLIFY = F)
  
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
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/3))
  }, b = b, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  Sgz2 <- mapply(function(b, S, SS, Fu, Fi, l0u, Delta, SigmaSplit){
    lhs.zeta <- Delta * S
    mu <- lapply(1:K, function(k) Fu[,b.inds[[k]],drop=F] %*% (gamma[k] * b[b.inds[[k]]]))
    mu <- Reduce('+', mu)
    temp <- sapply(1:K, function(k) crossprod(b[b.inds[[k]]]))
    Var <- lapply(1:K, function(k) sqrt(diag(gamma[k]^2 * tcrossprod(Fu[,b.inds[[k]]] %*% SigmaSplit[[k]], Fu[,b.inds[[k]]]))))
    tau <- Reduce('+', lapply(1:K, function(k) gamma[k]^2 * tcrossprod(Fu[,b.inds[[k]]] %*% SigmaSplit[[k]], Fu[,b.inds[[k]]])))
    tau <- sqrt(diag(tau))
    rhs.gamma <- matrix(0, nrow = gh, ncol = length(gamma))
    rhs.zeta <- matrix(0, nrow = gh, ncol = length(zeta))
    lhs.gamma <- Delta * sapply(1:K, function(k) Fi[,b.inds[[k]],drop=F] %*% b[b.inds[[k]]]); 
    for(g in 1:length(gamma)){
      for(l in 1:gh){
        rhs.gamma[l,g] <- w[l] * crossprod(l0u * exp(SS %*% zeta + gamma[g] * Fu %*% b + Var[[g]] * v[l]), Fu %*% b)
      }
    }
    
    for(l in 1:gh){
      rhs.zeta[l,] <- w[l] * crossprod(SS, l0u * exp(SS %*% zeta + Fu %*% (gamma.rep * b) + tau * v[l]))
    }
    
    
    c(lhs.gamma - c(colSums(rhs.gamma)), lhs.zeta - c(colSums(rhs.zeta)))
  }, b = b, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SigmaSplit = SigmaSplit,
     SIMPLIFY = F)
  
  Sg <- function()
  
  
  # Collate and form information --------------------------------------------
  Ss2 <- lapply(1:n, function(i){
    do.call(c, lapply(disps, function(j){
      Ss[[j]][[i]]
    }))
  })
  
  S <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, Sb, Ss, c(Sgz))
  }, sD = test2, Sb = Sb, Ss = Ss2, Sgz = Sgz)
  
  SS <- rowMeans(S)
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[,i]))) - SS
  
  sqrt(diag(solve(I)))
  
}