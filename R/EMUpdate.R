#' @keywords internal
EMupdate <- function(Omega, family, X, Y, Z, b,                # Longit.
                     S, SS, Fi, Fu, l0i, l0u, Delta, l0, sv,   # Survival
                     w, v, n, m, hessian,                      # Quadrature + additional info.
                     beta.inds, b.inds, K, q){
  
  #' Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
  beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) # Indexed for C++ use.
  gamma.rep <- rep(gamma, sapply(b.inds, length))
  
  #' Find b.hat and Sigma =====================
  if(hessian == 'auto') .hess <- T else .hess <- F
  b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
    optim(b, joint_density, joint_density_ddb,
          Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
          beta_inds = beta.inds2, b_inds = b.inds2, K = K,
          method = 'BFGS', hessian = .hess)
  }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
  Fu = Fu, l0u = l0u, SIMPLIFY = F)
  
  b.hat <- lapply(b.update, el, 1)
  if(.hess){
    Sigma <- lapply(b.update, function(x) solve(x$hessian))
  }else{
    Sigma <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
      solve(
        joint_density_sdb(b = b, Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family,
                          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep,
                          zeta = zeta, beta_inds = beta.inds2, b_inds = b.inds2, K = K, eps = 1e-5))
    }, b = b.hat, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
    Fu = Fu, l0u = l0u, SIMPLIFY = F)
  }
  
  # Split into K constituent sub-matrices.
  SigmaSplit <- lapply(Sigma, function(x) lapply(b.inds, function(y) as.matrix(x[y,y])))
  bsplit <- lapply(b.hat, function(x) lapply(b.inds, function(y) x[y]))                # Needed for updates to beta.
  bmat <- lapply(bsplit, bind.bs) # Needed for E[\ell(\gamma,\zeta)|\b...|\Omega].
  
  # #########
  # E-step ##
  # #########
  
  # D =========================================
  D.update <- mapply(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat, SIMPLIFY = F)
  
  # \beta =====================================
  Sb <- mapply(function(X, Y, Z, b){
    Sbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K)
  }, X = X, Y = Y, Z = Z, b = bsplit, SIMPLIFY = F)
  Hb <- mapply(function(X, Y, Z, b){
    Hbeta(beta, X, Y, Z, b, sigma, family, beta.inds2, K)
  }, X = X, Y = Y, Z = Z, b = bsplit, SIMPLIFY = F)
  
  # Dispersion ('\sigma') =====================
  funlist <- unlist(family)
  disps <- which(funlist %in% c('gaussian', 'genpois', 'Gamma'))
  sigma.update <- replicate(K, list(), simplify = F)
  for(j in disps){
    tau <- mapply(function(Z, S) unname(sqrt(diag(tcrossprod(Z[[j]] %*% S[[j]], Z[[j]])))), Z = Z, S = SigmaSplit)
    if(funlist[j] == 'gaussian'){
      sigma.update[[j]] <- sum(mapply(function(X, Y, Z, b, tau){
        vare_update(X[[j]], Y[[j]], Z[[j]], b[[j]], beta[beta.inds[[j]]], tau, w, v)
      }, X = X, Y = Y, Z = Z, b = bsplit, tau = tau))
    }
    if(funlist[j] == 'genpois'){
      sigma.update[[j]] <- rowSums(mapply(function(b, X, Y, Z, tau){
        unlist(phi_update(b[[j]], X[[j]], Y[[j]], Z[[j]], beta[beta.inds[[j]]], sigma[[j]], w, v, tau))
      }, b = bsplit, X = X, Y = Y, Z = Z, tau = tau))
    }
    if(funlist[j] == 'Gamma'){
      sigma.update[[j]] <- rowSums(mapply(function(X, Y, Z, tau, b){
        unlist(shape_update(sigma[[j]], X[[j]], Y[[j]], Z[[j]], tau, beta[beta.inds[[j]]], b[[j]], w, v))
      }, X = X, Y = Y, Z = Z, tau = tau, b = bsplit))
    }
  }
  for(j in setdiff(1:K, c(disps))) sigma.update[[j]] <- NA # Return null for all distsn which do not have disp. parameters.
  
  # Survival parameters (\gamma, \zeta) =======
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/3))
  }, b = b.hat, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta)
  
  Hgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, q, .Machine$double.eps^(1/4))
  }, b = b.hat, Sigma = SigmaSplit, S = S, SS = SS, Fu = Fu, Fi = Fi, l0u = l0u, Delta = Delta, SIMPLIFY = F)
  
  # #########
  # M-step ##
  # #########
  
  # D
  D.new <- Reduce('+', D.update)/n
  
  # beta
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  
  # Dispersion
  ms <- colSums(do.call(rbind, m))
  sigma.new <- as.list(rep(0., K))
  for(j in disps){ # This a little inefficient; kept as separating E- and M- steps
    if(funlist[j] == 'gaussian'){
      sigma.new[[j]] <- setNames(sigma.update[[j]]/ms[j],names(sigma[[j]]))
    }else if(funlist[j] == "genpois"){
      sigma.new[[j]] <- sigma[[j]] - sigma.update[[j]]['Score']/sigma.update[[j]]['Hessian']
    }else if(funlist[j] == 'Gamma'){
      sigma.new[[j]] <- sigma[[j]] - sigma.update[[j]]['Score']/sigma.update[[j]]['Hessian']
    }
  }
  
  # Survival parameters (gamma, zeta)
  gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz))
  
  # The baseline hazard and related objects
  lambda.update <- lambdaUpdate(sv$surv.times, do.call(cbind, sv$ft.mat), gamma, gamma.rep, zeta, S, SigmaSplit, b.hat, w, v, b.inds2, K, q)
  l0.new <- sv$nev/rowSums(lambda.update)
  l0u.new <- lapply(l0u, function(ll){
    l0.new[1:length(ll)]
  })
  l0i.new <- l0.new[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  
  # Return
  list(
    D = D.new, beta = beta.new, sigma = sigma.new,       # Yk responses
    gamma = gammazeta.new[1:K], zeta = gammazeta.new[(K+1):length(gammazeta.new)],  # Survival
    l0 = l0.new, l0u = l0u.new, l0i = as.list(l0i.new),  #   Hazard
    b = b.hat                                            #   REs.
  )
  
}