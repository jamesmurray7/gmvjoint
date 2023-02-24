.jd <- function(params, b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
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
                Fu, haz, gamma.rep, zeta, beta.inds2, b.inds2, K) * -1 - mvtnorm::dmvnorm(b, sigma = D, log = T) + 
    (-0.5 *log(2*pi) -0.5 * log(det(D)) - 0.5 * sum(diag(solve(D) %*% (Sigma + tcrossprod(b)))))
}

.jd(params,b = b[[1]], Y=Y[[1]], X=X[[1]],Z=Z[[1]],
    Delta = Delta[[1]], S=S[[1]], Fi=Fi[[1]],l0i=sv.new$l0i[[1]],SS = SS[[1]],
    Fu = Fu[[1]], haz = sv.new$l0u[[1]], Sigma=Sigma[[1]])

Scores <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
  pracma::grad(.jd, params, b = b, Y=Y, X=X,Z=Z,
               Delta = Delta, S=S, Fi=Fi,l0i=l0i,SS = SS,
               Fu = Fu, haz = haz, Sigma = Sigma)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = sv.new$l0i, 
SS = SS, Fu = Fu, haz = sv.new$l0u, Sigma = Sigma)

Hessian <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
  pracma::hessian(.jd, params, b = b, Y=Y, X=X,Z=Z,
                  Delta = Delta, S=S, Fi=Fi,l0i=l0i,SS = SS,
                  Fu = Fu, haz = haz, Sigma = Sigma)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = sv.new$l0i, 
SS = SS, Fu = Fu, haz = sv.new$l0u, Sigma = Sigma, SIMPLIFY = F)

.jd2 <- function(params, b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma,
                 C = TRUE, indiv = FALSE){
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
  
  if(!indiv){
    tau.long <- mapply(function(z, S) sqrt(diag(tcrossprod(z[[1]] %*% S, z[[1]]))), z = Z, S = Sigma, SIMPLIFY = F)
    tau.surv <- mapply(function(fu, S) sqrt(diag(gamma^2 * tcrossprod(fu %*% S, fu))), fu = Fu, S = Sigma, SIMPLIFY = F)
  }else{
    tau.long <- sqrt(diag(tcrossprod(Z[[1]] %*% Sigma, Z[[1]])))
    tau.surv <- sqrt(diag(gamma^2 * tcrossprod(Fu %*% S, Fu)))
  }
  if(!C){
    w <- 1;v <- 0;gh <- 1
    if(!indiv){
      Sigma <- lapply(Sigma, function(x) x * 0)
    }else{
      Sigma <- Sigma * 0
      tau.long <- tau.long * 0 
      tau.surv <- tau.surv * 0
    } 
  }
  
  if(!indiv){
    f.b <- mapply(function(b, S){
      -sv$q * 0.5 * log(2*pi) - 0.5 * log(det(D)) - 0.5 * sum(diag(solve(D) %*% (S + tcrossprod(b))))
    },b = b, S = Sigma)
    
    f.yb <- mapply(function(y, x, z, b, tau){
      eta <- x[[1]] %*% betas + z[[1]] %*% b
      f.yb <- numeric(length(y))
      for(l in 1:gh)
        f.yb <- f.yb + w[l] * dbinom(y[[1]], 1, plogis(eta + tau * v[l]), T)
      sum(f.yb)
    }, y = Y, x = X, z = Z, b = b, tau = tau.long)
    
    f.tb <- mapply(function(b, S, SS, fu, fi, l0u, l0i, tau, Del){
      lhs <- if(Del == 1) log(l0i) + S %*% zeta + fi %*% (gamma * b) else 0
      rhs <- numeric(gh)
      for(l in 1:gh)
        rhs[l] <- w[l] * crossprod(l0u, exp(SS %*% zeta + fu %*% (gamma * b) + tau * v[l]))
      lhs - sum(rhs)
    }, b = b, S = S, SS = SS, fu = Fu, fi = Fi, l0u = haz, l0i = l0i, tau =  tau.surv, Del = Delta)
  }else{
    f.b <- sv$q * 0.5 * log(2*pi) - 0.5 * log(det(D)) - 0.5 * sum(diag(solve(D) %*% (Sigma + tcrossprod(b))))
    eta <- X[[1]] %*% betas + Z[[1]] %*% b
    f.yb <- numeric(length(Y[[1]]))
    for(l in 1:gh)
      f.yb <- f.yb + w[l] * dbinom(Y[[1]], 1, plogis(eta + tau.long * v[l]), T)
    f.yb <- sum(f.yb)
    lhs <- if(Delta == 1) log(l0i) + S %*% zeta + Fi %*% (gamma * b) else 0
    rhs <- numeric(gh)
    for(l in 1:gh)
      rhs[l] <- w[l] * crossprod(haz, exp(SS %*% zeta + Fu %*% (gamma * b) + tau.surv * v[l]))
    f.tb <- lhs - sum(rhs)
  }
  
  sum(f.b + f.yb + f.tb)
}

.jd2(params,b = b, Y=Y, X=X,Z=Z,
     Delta = Delta, S=S, Fi=Fi,l0i=sv.new$l0i,SS = SS,
     Fu = Fu, haz = sv.new$l0u, Sigma=Sigma, FALSE, FALSE)
.jd2(params,b = b[[1]], Y=Y[[1]], X=X[[1]],Z=Z[[1]],
      Delta = Delta[[1]], S=S[[1]], Fi=Fi[[1]],l0i=sv.new$l0i[[1]],SS = SS[[1]],
      Fu = Fu[[1]], haz = sv.new$l0u[[1]], Sigma=Sigma[[1]], FALSE, TRUE)

ScoreCi <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
  pracma::grad(.jd2, params, b = b, Y=Y, X=X,Z=Z,
               Delta = Delta, S=S, Fi=Fi,l0i=l0i,SS = SS,
               Fu = Fu, haz = haz, Sigma = Sigma, C = TRUE, indiv = TRUE)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = sv.new$l0i,
  SS = SS, Fu = Fu, haz = sv.new$l0u, Sigma = Sigma)

sqrt(diag(
     solve(Reduce('+', lapply(1:n, function(x) tcrossprod(ScoreCi[,x]))) - tcrossprod(rowSums(ScoreCi))/n)
))

dAdO <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
  pracma::hessian(.jd2, params, b = b, Y=Y, X=X,Z=Z,
               Delta = Delta, S=S, Fi=Fi,l0i=l0i,SS = SS,
               Fu = Fu, haz = haz, Sigma = Sigma, C = TRUE, indiv = TRUE)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = sv.new$l0i,
SS = SS, Fu = Fu, haz = sv.new$l0u, Sigma = Sigma, SIMPLIFY = F)

AO <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
  pracma::grad(.jd2, params, b = b, Y=Y, X=X,Z=Z,
               Delta = Delta, S=S, Fi=Fi,l0i=l0i,SS = SS,
               Fu = Fu, haz = haz, Sigma = Sigma, C = FALSE, indiv = TRUE)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = sv.new$l0i,
SS = SS, Fu = Fu, haz = sv.new$l0u, Sigma = Sigma, SIMPLIFY = F)

SO <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, haz, Sigma){
  pracma::grad(.jd2, params, b = b, Y=Y, X=X,Z=Z,
               Delta = Delta, S=S, Fi=Fi,l0i=l0i,SS = SS,
               Fu = Fu, haz = haz, Sigma = Sigma, C = TRUE, indiv = TRUE)
}, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = sv.new$l0i,
SS = SS, Fu = Fu, haz = sv.new$l0u, Sigma = Sigma, SIMPLIFY = F)


temp <- mapply(function(dd, AO, SO){
  dd + AO %*% t(AO - SO)
}, dd = dAdO, AO = AO, SO = SO, SIMPLIFY = F)
tt <- Reduce('+', temp)

sqrt(diag(solve(-tt) %*% Reduce('+', lapply(SO, tcrossprod)) %*% t(solve(-tt))))

dAdO <- pracma::hessian(.jd2, params, b = b, Y=Y, X=X,Z=Z,
                         Delta = Delta, S=S, Fi=Fi,l0i=sv.new$l0i,SS = SS,
                         Fu = Fu, haz = sv.new$l0u, Sigma=Sigma, C = TRUE)
SO  <- pracma::grad(.jd2, params, b = b, Y=Y, X=X,Z=Z,
                    Delta = Delta, S=S, Fi=Fi,l0i=sv.new$l0i,SS = SS,
                    Fu = Fu, haz = sv.new$l0u, Sigma=Sigma, C = TRUE)
AO  <- pracma::grad(.jd2, params, b = b, Y=Y, X=X,Z=Z,
                    Delta = Delta, S=S, Fi=Fi,l0i=sv.new$l0i,SS = SS,
                    Fu = Fu, haz = sv.new$l0u, Sigma=Sigma, C = FALSE)
sqrt(diag(solve(-1 * (dAdO + (AO %*% t(AO - SO))))))
