#' @keywords internal
obs.emp.I <- function(Omega, dmats, surv, sv,
                      Sigma, SigmaSplit, b, bsplit, 
                      l0u, w, v, n, family, K, q, beta.inds, b.inds){
  # Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  sigma <- Omega$sigma
  gamma <- c(Omega$gamma)
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
  
  # Scores ------------------------------------------------------------------
  # The RE covariance matrix, D
  # Dinv <- solve(D)
  # postmult <- diag(1, nrow = nrow(Dinv), ncol = ncol(Dinv))
  # postmult[postmult == 0] <- 2 # off-diagonals have twice the contribution!
  # sD <- mapply(function(b, S){
  #   vech(t(0.5 * (Dinv %*% (S + tcrossprod(b)) %*% Dinv) - 0.5 * Dinv)) * vech(postmult)
  # }, b = b, S = Sigma, SIMPLIFY = F)
  
  # NB this same as commented segment above, just done "properly".
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
      dDEbbT <- delta.D[[i]] %*% EbbT
      term <- 0.5 * Dinv - 0.5 * Dinv %*% dDEbbT %*% Dinv
      out <- 0.5 * (t(-term) - term)
      sum(diag(out))
    },
    b = b, S = Sigma,
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
  Sgz <- mapply(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta){
    Sgammazeta_cd(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, b.inds2, K, .Machine$double.eps^(1/3))
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
  
  SS <- rowSums(S) # sum S
  #  observed empirical information matrix (Mclachlan and Krishnan, 2008).
  SiSiT <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i])))
  H <- SiSiT - tcrossprod(SS)/n
  
  return(list(Score = S,
              Hessian = H))
}

#' @keywords internal
obs.emp.I2 <- function(Omega, dmats, surv, sv, b, l0i, l0u,
                       w, v, n, family, K, q, beta.inds, b.inds){
  # Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  sigma <- Omega$sigma
  gamma <- c(Omega$gamma); gamma.rep <- rep(gamma, sapply(b.inds, length))
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
      dDEbbT <- delta.D[[i]] %*% EbbT
      term <- 0.5 * Dinv - 0.5 * Dinv %*% dDEbbT %*% Dinv
      out <- 0.5 * (t(-term) - term)
      sum(diag(out))
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
  
  SS <- rowSums(S) # sum S
  #  observed empirical information matrix (Mclachlan and Krishnan, 2008).
  SiSiT <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i])))
  H <- SiSiT - tcrossprod(SS)/n
  
  return(list(Score = S,
              Hessian = H,
              b.hat = b.hat,
              Sigma = Sigma))
}

#' Extract the variance-covariance matrix from a \code{joint} fit.
#' 
#' @details Uses the observed-empirical \strong{approximation} of information matrix 
#' (Mclachlan & Krishnan, 2008). The standard errors for the baseline hazard are not estimated. 
#' 
#' @param object a joint model fit by the \code{joint} function.
#' @param corr should the correlation matrix be returned instead of the variance-covariance?
#' @param ... extra arguments (none used).
#' 
#' @return A variance-covariance matrix for the joint model object.
#'
#' @author James Murray \email{j.murray7@@ncl.ac.uk}
#' 
#' @section Methodology: 
#' 
#' Many competing ways exist for obtaining the observed information matrix in an EM algorithm. 
#' In the context of joint modelling, the observed empirical approximation of the information 
#' matrix has been used previously (\code{joineRML}, Hickey et al. 2018). Elsewhere,
#' estimation of the observed information in a semi-parametric setting is outlined neatly in
#' Xu et al. (2014). Here, they advocate for approximation of this information matrix by 
#' numerical differentiation of the profile Fisher Score vector. We do not consider this 
#' methodology owing to its computational expense. That is, for each element of \eqn{\Omega} 
#' which is perturbed by some small amount \eqn{\tilde{\Omega}^{p}}, we must re-calculate
#' \eqn{\hat{b}_i} and \eqn{\hat{\Sigma}_i}.
#' 
#' @references 
#' 
#' Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. \code{joineRML}: a joint model and
#' software package for time-to-event and multivariate longitudinal outcomes.
#' \emph{BMC Med. Res. Methodol.} 2018; \strong{50}
#' 
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions.} Second Edition. 
#' Wiley-Interscience; 2008.
#' 
#' Xu C, Baines PD, Wang J. Standard error estimation using the EM algorithm for the joint 
#' modeling of survival and longitudinal data. \emph{Biostatistics} 2014; \strong{15}(4).
#' 
#' @method vcov joint
#' @export
#' 
#' @examples
#' # Univariate fit on PBC data -------------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify univariate fit
#' long.formulas <- list(
#'   albumin ~ time + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, family = list('gaussian'))
#' 
#' vcov(fit)
vcov.joint <- function(object, corr = FALSE, ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  
  v <- object$vcov
  if(corr) 
    return(cov2cor(v))
  else
    return(v)
}

