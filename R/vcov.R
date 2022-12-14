#' @keywords internal
obs.emp.I <- function(Omega, dmats, surv, sv,
                      Sigma, SigmaSplit, b, bsplit, 
                      l0u, w, v, n, family, K, q, beta.inds, b.inds, beta.quad){
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
      out <- 0.5 * (S + tcrossprod(b)) %*% (Dinv %*% delta.D[[i]] %*% Dinv)   
      lhs[i] + sum(diag(out))
    },
    b = b, S = Sigma,
    SIMPLIFY = T)
  }
  
  sD <- sapply(1:nrow(vech.indices), sDi)
  sD <- lapply(1:nrow(sD), function(x) sD[x, ]) # Cast to list
  
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
  I <- Reduce('+', lapply(1:n, function(i) tcrossprod(S[, i]))) - tcrossprod(SS)/n
  
  I
}

#' Extract the variance-covariance matrix from a \code{joint} fit.
#' 
#' @details Uses the observed-empirical \strong{approximation} of information matrix 
#' (Mclachlan & Krishnan, 2008). The estimates for the baseline hazard are not estimated. 
#' 
#' @param object a joint model fit by the \code{joint} function.
#' @param corr should the correlation matrix be returned instead of the variance-covariance?
#' @param ... extra arguments (none used).
#' 
#' @return A variance-covariance matrix for the joint model object.
#'
#' @author James Murray \email{j.murray7@@ncl.ac.uk}
#' @references 
#' 
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions.} Second Edition. 
#' Wiley-Interscience; 2008.
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

