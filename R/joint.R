#' Fit a joint model to time-to-event and multivariate longitudinal data
#'
#' @param long.formulas A list of formula objects specifying the \eqn{K} responses. Each must be 
#'        usable by \code{\link[glmmTMB]{glmmTMB}}. A restriction is that unique identifiers must 
#'        be named `id`, and increment in intervals of at exactly one.
#' @param surv.formula A formula specifying the time-to-event sub-model. Must be usable by 
#'   \code{\link[survival]{coxph}}.
#' @param data A \code{data.frame} containing all covariates and responses.
#' @param family A list of families corresponding in order to \code{long.formula}.
#' @param post.process Logical, should post processing be done to obtain standard errors and 
#'   log-likelihood? Defaults to \code{TRUE}.
#' @param control A list of control values: \describe{
#' 
#'   \item{\code{verbose}}{Logical: If \code{TRUE}, at each iteration parameter information will 
#'   be printed to console. Default is \code{verbose=FALSE}.}
#'   \item{\code{conv}}{Character: Either \code{"absolute"} or \code{"relative"} to invoke 
#'   absolute or relative difference as the convergence criterion. Default is 
#'   \code{conv="relative"}.}
#'   \item{\code{tol}}{Numeric: Tolerance value used to assess convergence. 
#'   Default is \code{tol=1e-2}.}
#'   \item{\code{correlated}}{Logical: Should covariance parameters \strong{between} responses 
#'   be estimated and used in determination of model convergence? Default is 
#'   \code{correlated=TRUE}. A choice of \code{correlated=FALSE} is equivalent to imposing the 
#'   belief that deviations in longitudinal trajectories are not correlated across responses, but
#'    can \strong{greatly decrease} computation time.}
#'   \item{\code{gh.nodes}}{Integer: Number of weights and abscissae to use in gauss-hermite 
#'   quadrature. Defaults to \code{gh.nodes=3}, which is usually sufficient.}
#'   \item{\code{gh.sigma}}{Numeric: Standard deviation for gauss-hermite approximation of normal
#'   distribution. Defaults to \code{gh.sigma=1}. This should rarely (if ever) need altering.}
#'   \item{\code{hessian}}{Character: Determines if the variance-covariance matrix for 
#'   \eqn{\hat{b}_i}, \eqn{\hat{\Sigma}_i} should be calculated as part of the \code{optim} step
#'   in minimising the negative log-likelihood, or calculated post-hoc using forward differencing.
#'   Default is \code{hessian="auto"} for the former, with \code{hessian="manual"} the 
#'   option for the latter.}
#'   \item{\code{return.inits}}{Logical: Should lists containing the initial conditions for the 
#'   longitudinal and survival sub-models be returned? Defaults to \code{return.inits=FALSE}.}
#'   \item{\code{beta.quad}}{Logical: Should gauss-hermite quadrature be used to appraise
#'   calculation of score and Hessian in updates to fixed effects \eqn{\beta}? Default is 
#'   \code{beta.quad=FALSE} which works very well in most situations. Dispersion parameters
#'   and survival pair are always calculated with quadrature.}
#' 
#' }
#' 
#' @returns An object with class \code{joint}. See \code{\link{joint.object}} for information. 
#' 
#' @details Function \code{joint} fits a joint model to time-to-event data and multivariate 
#' longitudinal data. The longitudinal data can be specified by numerous models encompassing
#' a fairly wide range of data. This joint model fit is achieved by the use of an approximate
#' EM algorithm first proposed in Bernhardt et al. (2015), and later used in the 'classic' 
#' multivariate joint model in Murray and Philipson (2022). Each longitudinal response is 
#' modelled by 
#' 
#' \deqn{h(E[Y_{ik}|b_{ik};\Omega]) = X_{ik}\beta_k + Z_{ik}b_{ik}} 
#' 
#' where \eqn{h} is a known, monotonic link function. An association is induced between the 
#' \eqn{K}th response and the hazard \eqn{\lambda_i(t)} by: 
#' 
#' \deqn{\lambda_i(t)=\lambda_0(t)\exp\{S_i^T\zeta + \sum_{k=1}^K\gamma_kW_k(t)b_{ik}\}} 
#' 
#' where \eqn{\gamma_k} is the association parameter and \eqn{W_k(t)} is the vector function of 
#' time imposed on the \eqn{K}th random effects structure (i.e. intercept-and-slope; spline).
#' 
#' @section Family specification:
#'   Currently, five families are available for implementation, spanning continuous, binary and 
#'   count data types: \describe{
#'   
#'     \item{\code{'gaussian'}}{Normally distributed. The identity link is used. A term 
#'     \eqn{\sigma_k} will be estimated, denoting the \emph{variance} of this response}
#'     \item{\code{'binomial'}}{For binary data types, a logit link is used.}
#'     \item{\code{'poisson'}}{For count data types where dispersion is either non-consequential 
#'     or ignored. A log link is used.}
#'     \item{\code{'genpois'}}{For count data types where dispersion is at least of some
#'     secondary interest. A log link is used. A term \eqn{\sigma_k} is estimated, denoting
#'     the dispersion, \eqn{\varphi} of the response. This follows interpretation of Zamani & 
#'     Ismail (2012): \eqn{\varphi>0}: Over-dispersion; \eqn{\varphi<0}: Under-dispersion.
#'     \eqn{Var[Y]=(1+\varphi)^2\mu}.}
#'     \item{\code{'Gamma'}}{For continuous data where a Gamma distribution might be sensible.
#'     The log link is used. A term \eqn{\sigma_k} is be estimated, denoting the shape of the
#'     distribution.}
#'   
#'   } 
#'   
#'   For families where dispersion is estimated, this is \strong{always} specified by an 
#'   "intercept-only" formula only. This might change in future.
#'   
#' @section Standard error estimation: 
#'   We follow the approximation of the observed empirical information matrix detailed by 
#'   Mclachlan and Krishnan (2008), and later used in \code{joineRML} (Hickey et al., 2018).
#'   These are only calculated if \code{post.process=TRUE}. Generally, these SEs are well-behaved,
#'   but their reliability will depend on multiple factors: Sample size; number of events; 
#'   collinearity of REs of responses; number of observed times, and so on.
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#'
#' @references 
#' 
#' Bernhardt PW, Zhang D and Wang HJ. A fast EM Algorithm for Fitting Joint Models of a Binary 
#' Response to Multiple Longitudinal Covariates Subject to Detection Limits. 
#' \emph{Computational Statistics and Data Analysis} 2015; \strong{85}; 37--53
#' 
#' Hickey GL, Philipson P, Jorgensen A, Kolamunnage-Dona R. \code{joineRML}: a joint model and
#' software package for time-to-event and multivariate longitudinal outcomes.
#' \emph{BMC Med. Res. Methodol.} 2018; \strong{50}
#'  
#' McLachlan GJ, Krishnan T. \emph{The EM Algorithm and Extensions.} Second Edition. 
#' Wiley-Interscience; 2008.
#' 
#' Murray, J and Philipson P. A fast approximate EM algorithm for joint models of survival and
#' multivariate longitudinal data.\emph{Computational Statistics and Data Analysis} 2022; 
#' \strong{170}; 107438
#' 
#' Zamani H and Ismail N. Functional Form for the Generalized Poisson Regression Model, 
#' \emph{Communications in Statistics - Theory and Methods} 2012; \strong{41(20)}; 3666-3675.
#' 
#' @seealso \code{\link{summary.joint}}, \code{\link{print.joint}}, \code{\link{fixef.joint}},
#' \code{\link{ranef.joint}}, \code{\link{vcov.joint}} and \code{\link{joint.object}}.
#' 
#' @export
#'
#' @examples
#' 
#' # 1) Fit on simulated bivariate data, (1x gaussian, 1x poisson) --------
#' beta <- do.call(rbind, replicate(2, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
#' gamma <- c(0.3, -0.3)
#' D <- diag(c(0.25, 0.09, 0.25, 0.05))
#' family <- list('gaussian', 'poisson')
#' data <- simData(ntms = 10, beta = beta, D = D, n = 100,
#'                 family = family, zeta = c(0, -0.2),
#'                 sigma = c(0.16, 0, 0.2), gamma = gamma)$data
#'
#' # Specify formulae and target families
#' long.formulas <- list(
#'   Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
#'   Y.2 ~ time + cont + bin + (1 + time|id)  # Poisson
#' )
#' surv.formula <- Surv(survtime, status) ~ bin
#' 
#' fit <- joint(long.formulas, surv.formula, data, family)
#' 
#' \donttest{
#' # 2) Fit on PBC data -----------------------------------------------------
#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'serBilir', 'albumin', 'spiders', 'platelets'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify GLMM sub-models, including interaction and natural spline terms
#' long.formulas <- list(
#'   serBilir ~ drug * (splines::ns(time, df = 3)) + (1 + splines::ns(time, df = 3)|id),
#'   albumin ~ drug * time + (1 + time|id),
#'   platelets ~ drug * time + (1 + time|id),
#'   spiders ~ drug * time + (1|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <-  joint(long.formulas, surv.formula, PBC, 
#'               family = list("gaussian", "gaussian", "poisson", "binomial"),
#'               control = list(verbose = T))
#' fit
#' }
joint <- function(long.formulas, surv.formula, data, family, post.process = TRUE, control = list()){
  
  start.time <- proc.time()[3]
  
  # Initial parsing ----
  formulas <- lapply(long.formulas, parseFormula)
  surv <- parseCoxph(surv.formula, data)
  n <- nrow(surv$survdata); K <- length(family)
  if(K!=length(long.formulas)) stop('Mismatched lengths of "family" and "long.formulas".')
  
  # Initial conditons ----
  inits.long <- Longit.inits(long.formulas, data, family)
  # Suss out indices of b_k and beta_k.
  b.inds <- lapply(1:K, function(k){
    nm <- inits.long$responses[k]
    which(grepl(nm, colnames(inits.long$b)))
  })
  beta.inds <- lapply(1:K, function(k){
    nm <- inits.long$responses[k]
    which(grepl(nm, names(inits.long$beta.init)))
  })
  q <- length(do.call(c, b.inds))
  
  inits.surv <- TimeVarCox(data, inits.long$b, surv$ph, formulas, b.inds)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  D <- inits.long$D.init
  sigma <- inits.long$sigma.init # dispersion / resid. variance / 0 otherwise.
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  # Survival parameters
  zeta <- inits.surv$inits[match(colnames(surv$ph$x), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma\\_', names(inits.surv$inits))]
  
  # Longitudinal and survival data objects ----
  dmats <- createDataMatrices(data, formulas)
  sv <- surv.mod(surv$ph, surv$survdata, formulas, inits.surv$l0.init)
  
  X <- dmats$X; Y <- dmats$Y; Z <- dmats$Z # Longitudinal data matrices
  m <- lapply(Y, function(y) sapply(y, length))
  # survival
  Fi <- sv$Fi; Fu <- sv$Fu; l0i <- sv$l0i; l0u <- sv$l0u; Delta <- surv$Delta 
  l0 <- sv$l0
  S <- sv$S; SS <- sv$SS
  
  # Assign family to joint density and parameter updates ----
  # Ensure family is a string, not a function i.e. if user supplies `Gamma` instead of `"Gamma"`.
  family <- lapply(family, function(family) if("function"%in%class(family)) family()$family else family)
  
  # Do we want to add-in correlated random effects between responses? For large K this greatly increases
  #   computation time and model instability.
  if(!is.null(control$correlated)) correlated <- control$correlated else correlated <- T
  if(!correlated) D[inits.long$off.inds] <- 0
  
  # Parameter vector and list ----
  Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, unlist(sigma)[inits.long$sigma.include], gamma, zeta)

  # Gauss-Hermite Quadrature ----
  if(!is.null(control$gh.nodes)) gh <- control$gh.nodes else gh <- 3
  if(!is.null(control$gh.sigma)) .sigma <- control$gh.sigma else .sigma <- 1
  
  GH <- statmod::gauss.quad.prob(gh, 'normal', sigma = .sigma)
  w <- GH$w; v <- GH$n
  
  # Begin EM ----
  diff <- 100; iter <- 0;
  if(!is.null(control$tol)) tol <- control$tol else tol <- 1e-2
  if(!is.null(control$maxit)) maxit <- control$maxit else maxit <- 200
  if(!is.null(control$conv)) conv <- control$conv else conv <- "relative"
  if(!conv%in%c('absolute', 'relative')) stop('Only "absolute" and "relative" difference convergence criteria implemented.')
  if(!is.null(control$verbose)) verbose <- control$verbose else verbose <- F
  if(!is.null(control$hessian)) hessian <- control$hessian else hessian <- 'auto'
  if(!hessian %in% c('auto', 'manual')) stop("Argument 'hessian' needs to be either 'auto' (i.e. from optim) or 'manual' (i.e. from _sdb, the defualt).")
  if(!is.null(control$return.inits)) return.inits <- control$return.inits else return.inits <- F
  if(!is.null(control$beta.quad)) beta.quad <- control$beta.quad else beta.quad <- F
  
  EMstart <- proc.time()[3]
  while(diff > tol && iter < maxit){
    update <- EMupdate(Omega, family, X, Y, Z, b, 
                       S, SS, Fi, Fu, l0i, l0u, Delta, l0, sv, 
                       w, v, n, m, hessian, beta.inds, b.inds, K, q, beta.quad)
    if(!correlated) update$D[inits.long$off.inds] <- 0
    params.new <- c(vech(update$D), update$beta, unlist(update$sigma)[inits.long$sigma.include], 
                    update$gamma, update$zeta)
    names(params.new) <- names(params)
    
    # Check (& print) convergence.
    diffs <- difference(params, params.new, conv)
    diff <- max(diffs)
    if(verbose){
      print(sapply(params.new, round, 4))
      message("Iteration ", iter + 1, ' maximum ', conv, ' difference: ', round(diff, 4), ' for parameter ', names(params.new)[which.max(diffs)])
      message("Largest RE relative difference: ", round(max(difference(do.call(cbind, update$b), do.call(cbind, b), conv)), 4))
    }
    
    # Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; sigma <- update$sigma
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    iter <- iter + 1
    Omega <- list(D = D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
    params <- params.new
  }
  
  EMend <- proc.time()[3]
  coeffs <- Omega
  coeffs$beta <- setNames(c(Omega$beta), names(inits.long$beta.init))
  out <- list(coeffs = coeffs,
              hazard = cbind(ft = sv$ft, haz = l0, nev = sv$nev))
  
  ModelInfo <- list()
  ModelInfo$ResponseInfo <- sapply(1:K, function(k){
    paste0(inits.long$responses[k], ' (', family[k], ')')
  })
  ModelInfo$family <- family
  ModelInfo$long.formulas <- long.formulas
  ModelInfo$surv.formulas <- surv.formula
  ModelInfo$control <- if(!is.null(control)) control else NULL
  ModelInfo$inds <- list(beta = beta.inds, b = b.inds)
  ModelInfo$nobs <- colSums(do.call(rbind, m))
  ModelInfo$n <- n
  ModelInfo$nev <- sum(sv$nev)
  out$ModelInfo <- ModelInfo
  
  # Post processing ----
  if(post.process){
    if(verbose) message('\nCalculating SEs...')
    beta.inds2 <- lapply(beta.inds, function(x) x - 1); b.inds2 <- lapply(b.inds, function(x) x - 1) 
    gamma.rep <- rep(gamma, sapply(b.inds, length))
    pp.start.time <- proc.time()[3]
    
    # b and Sigmai at MLEs
    b.update <- mapply(function(b, Y, X, Z, Delta, S, Fi, l0i, SS, Fu, l0u){
      optim(b, joint_density, joint_density_ddb,
            Y = Y, X = X, Z = Z, beta = beta, D = D, sigma = sigma, family = family, 
            Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS, Fu = Fu, haz = l0u, gamma_rep = gamma.rep, zeta = zeta,
            beta_inds = beta.inds2, b_inds = b.inds2, K = K,
            method = 'BFGS', hessian = T)
    }, b = b, Y = Y, X = X, Z = Z, Delta = Delta, S = S, Fi = Fi, l0i = l0i, SS = SS,
    Fu = Fu, l0u = l0u, SIMPLIFY = F)
    Sigma <- lapply(b.update, function(x) solve(x$hessian))
    b <- lapply(b.update, function(x) x$par)
    SigmaSplit <- lapply(Sigma, function(x) lapply(b.inds, function(y) as.matrix(x[y,y])))
    bsplit <- lapply(b, function(x) lapply(b.inds, function(y) x[y])) # Needed for updates to beta.
    # The Information matrix
    I <- structure(obs.emp.I(coeffs, dmats, surv, sv, 
                             Sigma, SigmaSplit, b, bsplit, 
                             l0u, w, v, n, family, K, q, beta.inds, b.inds, beta.quad),
                   dimnames = list(names(params), names(params)))

    I.inv <- tryCatch(solve(I), error = function(e) e)
    if(inherits(I.inv, 'error')) I.inv <- structure(MASS::ginv(I),
                                                    dimnames = dimnames(I))
    out$SE <- sqrt(diag(I.inv))
    out$vcov <- I.inv
    
    postprocess.time <- round(proc.time()[3]-pp.start.time, 2)
    # Calculate log-likelihood. Done separately as EMtime + postprocess.time is for EM + SEs.
    out$logLik <- joint.log.lik(coeffs, dmats, b, surv, sv, l0u, l0i, gamma.rep, beta.inds, b.inds, 
                                K, q, family, Sigma)
    # Collate RE and their variance
    REs <- do.call(rbind, b)
    attr(REs, 'Var') <- do.call(rbind, lapply(Sigma, diag))
    out$REs <- REs
  }
  comp.time <- round(proc.time()[3] - start.time, 3)
  out$elapsed.time <- c(`EM time` = unname(round(EMend - EMstart, 3)),
                        `Post processing` = if(post.process) unname(postprocess.time) else NULL,
                        `Total Computation time` = unname(comp.time),
                        `iterations` = iter)
  
  if(return.inits) out$inits = inits.long
  class(out) <- 'joint'
  return(out)
}

#' Print a \code{joint} object
#' @method print joint
#' @keywords internal
#' @export
print.joint <- function(x, ...){
  if(!inherits(x, 'joint')) stop('x must be a "joint" object!')
  
  M <- x$ModelInfo
  K <- length(M$ResponseInfo) # Number of responses
  
  # Data information
  cat(sprintf("Number of subjects: %d\n", M$n))
  cat(sprintf("Number of events: %d (%.2f%%)\n", M$nev, 100 * M$nev/M$n))
  
  # Longitudinal information: 
  cat("\n===================\nModel specification\n===================\n")
  if(K == 1){
    cat("Univariate longitudinal process specification:\n")
    cat(sprintf("%s: %s\n", M$ResponseInfo[1], deparse(M$long.formulas[[1]])))
  }else{
    cat("Multivariate longitudinal process specifications: \n")
    for(k in 1:K){
      cat(sprintf("%s: %s\n", M$ResponseInfo[k], deparse(M$long.formulas[[k]])))
    }
  }
  
  cat("\nSurvival sub-model specification: \n")
  cat(deparse(M$surv.formulas))
  
  cat("\n\nAssociation parameter estimates: \n")
  print(setNames(x$coeffs$gamma,
        M$ResponseInfo))
  cat("\n")
  invisible(1+1)
}
