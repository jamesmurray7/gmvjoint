#' @keywords internal
obs.emp.I <- function(Omega, dmats, surv, sv, family,
                      b, l0i, l0u, w, v, inds){
  # Unpack Omega ----
  D <- Omega$D
  beta <- c(Omega$beta)
  sigma <- Omega$sigma
  gamma <- c(Omega$gamma); gamma.rep <- rep(gamma, sapply(inds$R$b, length))
  zeta <- c(Omega$zeta)
  
  # Calculate b.hat, Sigma.hat at MLEs
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
  
  SigmaSplit <- lapply(Sigma, function(x) lapply(inds$R$b, function(y) as.matrix(x[y,y])))
  
  # Profile estimate for (gamma, zeta) at MLEs and b.hat
  l0.hat <- lambda_hat(b.hat, sv$Fu, sv$SS, Sigma, gamma.rep, zeta, sv$nev, w, v) # Profile estimate to work out (gamma, zeta).
  l0u.hat <- lapply(sv$l0u, function(ll){
    l0.hat[1:length(ll)]
  })
  
  # Scores ------------------------------------------------------------------
  # D =========================================
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
  
  # \beta =====================================
  tau <- Map(make_tau, Z = dmats$Z, S = SigmaSplit)
  Sb <- Map(function(X, Y, Z, W, b, tau){
    c(Sbeta(beta, X, Y, Z, W, b, sigma, family, inds$Cpp$beta, inds$Cpp$b, 
          dmats$K, tau, w, v) )
  }, X = dmats$X, Y = dmats$Y, Z = dmats$Z, W = dmats$W, b = b.hat, tau = tau)
  
  # Dispersion ('\sigma') =====================
  funlist <- unlist(family)
  eta <- Map(function(X, Z, b) make_eta(X, Z, beta, b, inds$Cpp$beta, inds$Cpp$b), 
             X = dmats$X, Z = dmats$Z, b = b.hat)
  disps <- which(funlist %in% c('gaussian', 'genpois', 'Gamma', 'negbin'))
  
  # Scores in form [[k]][[i]]
  Ssig <- lapply(seq_along(funlist), function(f){
    
    tau.f <- lapply(tau, el, f)
    eta.f <- lapply(eta, el, f)
    ff <- funlist[f]
    
    if(ff == "gaussian"){
      return(Map(function(Y, eta, tau, mi){
        rhs <- numeric(length(w))
        for(l in 1:length(w)){
          rhs[l] <- w[l] * crossprod(Y[[f]] - eta - tau * v[l])
        }
        -mi[[f]]/(2 * sigma[[f]]) + sum(rhs)/(2 * sigma[[f]]^2)
      }, Y = dmats$Y, eta = eta.f, tau = tau.f, mi = dmats$mi))
    }else if(ff == "Gamma"){
      return(Map(function(eta, Y, tau, W){
        pracma::grad(appxE_Gammasigma, sigma[[f]],
                     eta = eta, Y = Y[[f]], tau = tau, W = W[[f]], 
                     w = w, v = v)
      }, eta = eta.f, Y = dmats$Y, tau = tau.f, W = dmats$W))
    }else if(ff == "negbin"){
      return(Map(function(eta, Y, tau, W){
        pracma::grad(appxE_NegBinsigma, sigma[[f]],
                     eta = eta, Y = Y[[f]], tau = tau, W = W[[f]],
                     w = w, v = v)
      }, eta = eta.f, Y = dmats$Y, tau = tau.f, W = dmats$W))
    }else if(ff == "genpois"){
      return(Map(function(eta, Y, tau, W){
        pracma::grad(appxE_GenPoissigma, sigma[[f]],
                     eta = eta, Y = Y[[f]], tau = tau, W = W[[f]],
                     w = w, v = v)
      }, eta = eta.f, Y = dmats$Y, tau = tau.f, W = dmats$W))
    }else{
      return(NULL)
    }
    
  })
  
  # and swap around to obtain scores in form [[i]][[k]]
  Ssig <- lapply(1:dmats$n, function(i) lapply(1:dmats$K, function(k) Ssig[[k]][[i]]))
  
  # Survival parameters (\gamma, \zeta) =======
  Psurv <- length(c(gamma, zeta))
  Sgz <- Map(function(b, Sigma, S, SS, Fu, Fi, l0u, Delta, ST){
    if(length(ST))
      return(Sgammazeta(c(gamma, zeta), b, Sigma, S, SS, Fu, Fi, l0u, Delta, w, v, 
                        inds$Cpp$b, dmats$K, .Machine$double.eps^(1/3)))
    else
      return(rep(0, Psurv))
  }, b = b.hat, Sigma = Sigma, S = sv$S, SS = sv$SS, Fu = sv$Fu, 
  Fi = sv$Fi, l0u = l0u.hat, Delta = surv$Delta, ST = sv$surv.times)
  
  # Collate and form information --------------------------------------------
  S <- mapply(function(sD, Sb, Ss, Sgz){
    c(sD, Sb, unlist(Ss), Sgz)
  }, sD = sD, Sb = Sb, Ss = Ssig, Sgz = Sgz)
  
  SS <- rowSums(S) # sum 
  #  observed empirical information matrix (Mclachlan and Krishnan, 2008).
  SiSiT <- Reduce('+', lapply(1:dmats$n, function(i) tcrossprod(S[, i])))
  H <- SiSiT - tcrossprod(SS)/dmats$n
  
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

