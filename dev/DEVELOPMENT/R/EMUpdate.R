#' @keywords internal
EMupdate <- function(Omega, family, dmats, b,                # Params; families; dmats; REs;
                     sv, surv, w, v,                         # Survival; quadrature;
                     con, inds){                             # control arguments; indices
  
  # Unpack Omega, the parameter vector
  D <- Omega$D; beta <- Omega$beta; sigma <- Omega$sigma; gamma <- Omega$gamma; zeta <- Omega$zeta
  gamma.rep <- rep(gamma, sapply(inds$R$b, length))
  
  # Find b.hat and Sigma =====================
  b.update <- Map(function(b, Y, X, Z, W, Delta, S, Fi, l0i, SS, Fu, l0u){
    optim(b, joint_density, joint_density_ddb,
          Y = Y, X = X, Z = Z, W = W, beta = beta, D = D, sigma = sigma, family = family, 
          Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
          beta_inds = inds$Cpp$beta, b_inds = inds$Cpp$b, K = dmats$K,
          method = 'BFGS', hessian = TRUE)
  }, b = b, Y = dmats$Y, X = dmats$X, Z = dmats$Z, W = dmats$W, Delta = surv$Delta, 
     S = sv$S, Fi = sv$Fi, l0i = sv$l0i, SS = sv$SS, Fu = sv$Fu, l0u = sv$l0u)
  
  b.hat <- lapply(b.update, el, 1)
  Sigma <- lapply(b.update, function(x) solve(x$hessian))
  
  # Split into K constituent sub-matrices.
  SigmaSplit <- lapply(Sigma, function(x) lapply(inds$R$b, function(y) as.matrix(x[y,y])))
  
  # #########
  # E-step ##
  # #########
  
  # D =========================================
  D.update <- Map(function(Sigma, b) Sigma + tcrossprod(b), Sigma = Sigma, b = b.hat)
  
  # \beta =====================================
  tau <- Map(make_tau, Z = dmats$Z, S = SigmaSplit)
  eta <- Map(function(X, Z, b) make_eta(X, Z, beta, b, inds$Cpp$beta, inds$Cpp$b), 
             X = dmats$X, Z = dmats$Z, b = b.hat)
  
  Sb <- Map(function(X, Y, Z, W, b, tau){
    Sbeta(beta, X, Y, Z, W, b, sigma, family, inds$Cpp$beta, inds$Cpp$b, 
          dmats$K, tau, w, v) 
  }, X = dmats$X, Y = dmats$Y, Z = dmats$Z, W = dmats$W, b = b.hat, tau = tau)
  
  Hb <- Map(function(X, Y, Z, W, b, tau){
    Hbeta(beta, X, Y, Z, W, b, sigma, family, inds$Cpp$beta, inds$Cpp$b, 
          dmats$K, tau, w, v) 
  }, X = dmats$X, Y = dmats$Y, Z = dmats$Z, W = dmats$W, b = b.hat, tau = tau)
  
  # Dispersion ('\sigma') =====================
  funlist <- unlist(family)
  disps <- which(funlist %in% c('gaussian', 'genpois', 'Gamma', 'negbin'))
  
  sigma.update <- lapply(seq_along(funlist), function(f){
    
    tau.f <- lapply(tau, el, f)
    eta.f <- lapply(eta, el, f)
    ff <- funlist[f];
    
    if(ff == "gaussian"){
      return(
        sum(mapply(function(eta, Y, tau){
          sigma2_Gaussian_update(eta, Y[[f]], tau, w, v)
        }, eta = eta.f, Y = dmats$Y, tau = tau.f))/dmats$m[f]
      )
    }else if(ff == "Gamma"){
      S_H_ <- Map(function(eta, Y, tau, W){
        S <- pracma::grad(appxE_Gammasigma, sigma[[f]],
                          eta = eta, Y = Y[[f]], tau = tau, W = W[[f]], w = w, v = v)
        H <- pracma::hessian(appxE_Gammasigma, sigma[[f]],
                             eta = eta, Y = Y[[f]], tau = tau, W = W[[f]], w = w, v = v)
        list(S = S, H = H)
      }, eta = eta.f, Y = dmats$Y, tau = tau.f, W = dmats$W)
      return(solve(Reduce('+', lapply(S_H_, el, 2)), Reduce('+', lapply(S_H_, el, 1))))
    }else if(ff == "negbin"){
      S_H_ <- Map(function(eta, Y, tau, W){
        S <- pracma::grad(appxE_NegBinsigma, sigma[[f]],
                          eta = eta, Y = Y[[f]], tau = tau, W = W[[f]],
                          w = w, v = v)
        H <- pracma::hessian(appxE_NegBinsigma, sigma[[f]],
                             eta = eta, Y = Y[[f]], tau = tau, W = W[[f]],
                             w = w, v = v)
        list(S = S, H = H)
      }, eta = eta.f, Y = dmats$Y, tau = tau.f, W = dmats$W)
      return(solve(Reduce('+', lapply(S_H_, el, 2)), Reduce('+', lapply(S_H_, el, 1))))
    }else if(ff == "genpois"){
      S_H_ <- Map(function(eta, Y, tau, W){
        S <- pracma::grad(appxE_GenPoissigma, sigma[[f]],
                          eta = eta, Y = Y[[f]], tau = tau, W = W[[f]],
                          w = w, v = v)
        H <- pracma::hessian(appxE_GenPoissigma, sigma[[f]],
                             eta = eta, Y = Y[[f]], tau = tau, W = W[[f]],
                             w = w, v = v)
        list(S = S, H = H)
      }, eta = eta.f, Y = dmats$Y, tau = tau.f, W = dmats$W)
      return(solve(Reduce('+', lapply(S_H_, el, 2)), Reduce('+', lapply(S_H_, el, 1))))
    }else{
      return(0)
    }
    
  })
  
  # Survival parameters (\gamma, \zeta) =======
  l0.hat <- lambda_hat(b.hat, sv$Fu, sv$SS, Sigma, gamma.rep, zeta, sv$nev, w, v) # Profile estimate to work out (gamma, zeta).
  l0u.hat <- lapply(sv$l0u, function(ll){                                         # _Not_ estimator for lambda_0.
    l0.hat[1:length(ll)]
  })
  
  Psurv <- length(c(gamma, zeta))
  
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta, ST){
    if(length(ST))
      return(Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 
                        inds$Cpp$b, dmats$K, .Machine$double.eps^(1/3)))
    else
      return(rep(0, Psurv))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
     Fi = sv$Fi, l0u = l0u.hat, Delta = surv$Delta, ST = sv$surv.times)
  
  Hgz <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta, ST){
    if(length(ST))
      return(Hgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 
                          inds$Cpp$b, dmats$K, .Machine$double.eps^(1/4)))
    else
      return(matrix(0, Psurv, Psurv))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
     Fi = sv$Fi, l0u = l0u.hat, Delta = surv$Delta, ST = sv$surv.times)
  
  # #########
  # M-step ##
  # #########
  
  # D
  D.new <- Reduce('+', D.update)/dmats$n
  
  # beta
  beta.new <- beta - solve(Reduce('+', Hb), Reduce('+', Sb))
  
  # Dispersion
  sigma.new <- replicate(dmats$K, 0, simplify = FALSE)
  for(j in disps){ # This a little inefficient; kept as separating E- and M- steps
    if(funlist[j] == 'gaussian'){
      sigma.new[[j]] <- setNames(sigma.update[[j]], names(sigma[[j]]))
    }else{
      sigma.new[[j]] <- sigma[[j]] - sigma.update[[j]]
    }
  }

  # Survival parameters (gamma, zeta)
  gammazeta.new <- c(gamma, zeta) - solve(Reduce('+', Hgz), rowSums(Sgz))
  gamma.new <- gammazeta.new[1:dmats$K]; zeta.new <- gammazeta.new[(dmats$K+1):length(gammazeta.new)]
  
  # The baseline hazard and related objects
  lambda.update <- c(lambda_update(b.hat, sv$Fu, sv$SS, Sigma, sv$surv.times, 
                                   rep(gamma.new, sapply(inds$R$b, length)), zeta.new, sv$nev,
                                   w, v))
  
  l0u.new <- lapply(sv$l0u, function(ll){
    lambda.update[1:length(ll)]
  })
  l0i.new <- lambda.update[match(sv$Ti, sv$ft)] 
  l0i.new[is.na(l0i.new)] <- 0
  
  # Return
  list(
    D = D.new, beta = beta.new, sigma = sigma.new,              #   K responses;
    gamma = gamma.new, zeta = zeta.new,                         #   survival;
    l0 = lambda.update, l0u = l0u.new, l0i = as.list(l0i.new),  #   hazard;
    l0u.hat = l0u.hat,                                          #   (+ profile);
    b = b.hat,                                                  #   REs;
    Sigma = Sigma                                               #   + their variance.
  )
  
}