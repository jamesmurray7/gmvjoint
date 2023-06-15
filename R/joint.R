#' Fit a joint model to time-to-event and multivariate longitudinal data
#'
#' @param long.formulas A list of formula objects specifying the \eqn{K} responses. Each must be 
#'        usable by \code{\link[glmmTMB]{glmmTMB}}. A restriction is that unique identifiers must 
#'        be named \code{id}, and increment in intervals of at exactly one. The variable for time
#'        must be named \code{time}.
#' @param surv.formula A formula specifying the time-to-event sub-model. Must be usable by 
#'   \code{\link[survival]{coxph}}.
#' @param data A \code{data.frame} containing all covariates and responses.
#' @param family A list of length \eqn{K} containing strings denoting the exponential families 
#' for each longitudinal sub-model, corresponding in order to \code{long.formulas}. For choices 
#' of \code{family}, see \strong{details}.
#' @param disp.formulas An optional list of length \eqn{K} specifying the dispersion models
#' wanted for each longitudinal sub-model, corresponding in order to \code{long.formulas}. Defaults
#' to \code{disp.formulas = NULL}. See \strong{details} for more information.
#' @param control A list of control values: \describe{
#' 
#'   \item{\code{verbose}}{Logical: If \code{TRUE}, at each iteration parameter information will 
#'   be printed to console. Default is \code{verbose=FALSE}.}
#'   \item{\code{conv}}{Character: Convergence criterion, see \strong{details}.}
#'   \item{\code{tol.abs}}{Numeric: Tolerance value used to assess convergence, see 
#'   \strong{details}. Default is \code{tol.abs=1e-3}.}
#'   \item{\code{tol.rel}}{Numeric: Tolerance value used to assess convergence, see 
#'   \strong{details}. Default is \code{tol.rel=1e-2}.}
#'   \item{\code{tol.den}}{Numeric: Tolerance value used to assess convergence, see 
#'   \strong{details}. Default is \code{tol.den=1e-3}.}
#'   \item{\code{tol.thr}}{Numeric: Threshold used when \code{conv = 'sas'}, see 
#'   \strong{details}. Default is \code{tol.thr=1e-1}.}
#'   \item{\code{maxit}}{Integer: Maximum number of EM iterations to carry out before
#'   exiting the algorithm. Defaults to \code{maxit=200L}, which is usually sufficient.}
#'   \item{\code{correlated}}{Logical: Should covariance parameters \strong{between} responses 
#'   be estimated and used in determination of model convergence? Default is 
#'   \code{correlated=TRUE}. A choice of \code{correlated=FALSE} is equivalent to imposing the 
#'   belief that deviations in longitudinal trajectories are not correlated across responses, but
#'   can decrease computation time, particularly for large \eqn{K}.}
#'   \item{\code{gh.nodes}}{Integer: Number of weights and abscissae to use in gauss-hermite 
#'   quadrature. Defaults to \code{gh.nodes=3}, which is usually sufficient.}
#'   \item{\code{gh.sigma}}{Numeric: Standard deviation for gauss-hermite approximation of normal
#'   distribution. Defaults to \code{gh.sigma=1}. This should rarely (if ever) need altering.}
#'   \item{\code{return.dmats}}{Logical: Should data matrices be returned? Defaults to 
#'   \code{return.dmats=TRUE}. Note that some S3 methods for \code{\link{joint.object}}s
#'   require the returned object to include these data matrices.}
#'   \item{\code{return.inits}}{Logical: Should a list of inital conditons be returned? 
#'   Defaults to \code{return.inits=FALSE}.}
#'   \item{\code{center.ph}}{Logical: Should the survival covariates be mean-centered? Defaults
#'   to \code{center.ph=TRUE}.}
#'   \item{\code{post.process}}{Logical: Should model post-processing be carried out (assumes
#'   that the model has converged). Defaults to \code{post.process = TRUE} which then returns
#'   posterior modes and their variance for the random effects, as well as approximated standard
#'   error. This is largely for internal use (i.e. if bootstrapping to obtain SEs instead).}
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
#' \deqn{h_k(E[Y_{ik}|b_{ik};\Omega]) = X_{ik}\beta_k + Z_{ik}b_{ik}} 
#' 
#' where \eqn{h_k} is a known, monotonic link function. An association is induced between the 
#' \eqn{K}th response and the hazard \eqn{\lambda_i(t)} by: 
#' 
#' \deqn{\lambda_i(t)=\lambda_0(t)\exp\{S_i^T\zeta + \sum_{k=1}^K\gamma_kW_k(t)^Tb_{ik}\}} 
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
#'     The log link is used. A term \eqn{\sigma_k} is be estimated, denoting the (log) shape of 
#'     the distribution, which is then reported as \eqn{\varphi_k=\exp\{\sigma_k\}}.}
#'     \item{\code{"negbin"}}{For count data types where overdispersion is modelled. A log link
#'     is used. A term \eqn{\sigma_k} is estimated, which is then reported as 
#'     \eqn{\varphi_k=\exp\{\sigma_k\}} which is the overdispersion. The variance of the response
#'     is \eqn{Var[Y]=\mu+\mu^2/\varphi}.}
#'   
#'   } 
#'   
#'  For families \code{"negbin"}, \code{"Gamma"}, \code{"genpois"}, the user can define the 
#'  dispersion model desired in \code{disp.formulas}. For the \code{"negbin"} and \code{"Gamma"}
#'  cases, we define \eqn{\varphi_i=\exp\{W_i\sigma_i\}} (i.e. the exponent of the linear 
#'  predictor of the dispersion model; and for \code{"genpois"} the identity of the linear
#'  is used.
#' 
#' @section Dispersion models:
#'   The \code{disp.formulas} in the function call allows the user to model the dispersion for
#'   a given sub-model if wanted. The default value \code{disp.formulas = NULL} simply imposes
#'   an 'intercept only' model. If the \eqn{k}th item in \code{disp.formulas} corresponds to 
#'   a longitudinal sub-model with no dispersion term, then it is simply ignored. With this in 
#'   mind then, if a dispersion model is only required for, say, one sub-model, then the 
#'   corresponding item in this list of models should be specified as such, with the others set to
#'   \code{~1}.
#'
#' @section Standard error estimation: 
#'   We follow the approximation of the observed empirical information matrix detailed by 
#'   Mclachlan and Krishnan (2008), and later used in \code{joineRML} (Hickey et al., 2018).
#'   These are only calculated if \code{post.process=TRUE}. Generally, these SEs are well-behaved,
#'   but their reliability will depend on multiple factors: Sample size; number of events; 
#'   collinearity of REs of responses; number of observed times, and so on. Some more discussion/
#'   references are given in \code{\link{vcov.joint}}.
#'   
#' @section Convergence of the algorithm:
#' A few convergence criteria (specified by \code{control$conv}) are available: \describe{
#'   \item{\code{abs}}{Convergence reached when maximum absolute change in parameter estimates
#'   is \code{<tol.abs}.}
#'   \item{\code{rel}}{Convergence reached when maximum absolute relative change in parameter
#'   estimates is \code{<tol.rel}. A small amount (\code{tol.den}) is added to the denominator 
#'   to eschew numerical issues if parameters are nearly zero.}
#'   \item{\code{either}}{Convergence is reached when either \code{abs} or \code{rel} are met.}
#'   \item{\code{sas}}{Assess convergence for parameters \eqn{|\Omega_a|}\code{<tol.thr} by the
#'   \code{abs} criterion, else \code{rel}. This is the default.}
#' 
#' }
#' Note that the baseline hazard is updated at each EM iteration, but is not monitored for 
#' convergence.
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
#' @seealso \code{\link{summary.joint}}, \code{\link{logLik.joint}}, 
#' \code{\link{extractAIC.joint}}, \code{\link{fixef.joint}}, \code{\link{ranef.joint}},
#' \code{\link{vcov.joint}}, \code{\link{joint.object}} and \code{\link{xtable.joint}}. For
#' data simulation see \code{\link{simData}}.
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
#'                 sigma = list(0.16, 0), gamma = gamma)$data
#'
#' # Specify formulae and target families
#' long.formulas <- list(
#'   Y.1 ~ time + cont + bin + (1 + time|id),  # Gaussian
#'   Y.2 ~ time + cont + bin + (1 + time|id)   # Poisson
#' )
#' surv.formula <- Surv(survtime, status) ~ bin
#' 
#' fit <- joint(long.formulas, surv.formula, data, family)
#' 
#' \donttest{
#' # 2) Fit on PBC data -----------------------------------------------------
#' data(PBC)
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'serBilir', 'albumin', 'spiders', 'platelets'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify GLMM sub-models, including interaction and quadratic time terms
#' long.formulas <- list(
#'   log(serBilir) ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id),
#'   albumin ~ drug * time + (1 + time|id),
#'   platelets ~ drug * time + (1 + time|id),
#'   spiders ~ drug * time + (1|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <-  joint(long.formulas, surv.formula, PBC, 
#'               family = list("gaussian", "gaussian", "poisson", "binomial"),
#'               control = list(verbose = TRUE))
#' fit
#' }
#' \donttest{
#' # 3) Fit with dispersion models ----------------------------------------
#' beta <- do.call(rbind, replicate(2, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
#' gamma <- c(0.3, -0.3)
#' D <- diag(c(0.25, 0.09, 0.25, 0.05))
#' family <- list('negbin', 'poisson')   # As an example; only requires one dispersion model.
#' sigma <- list(c(1, 0.2), 0)           # Need to specify the model in simData call too.
#' disp.formulas = list(~time, ~1)       # Even though poisson doesn't model dispersion, need to
#'                                       # populate this element in disp.formulas!
#' # Simulate some data
#' data <- simData(ntms = 10, beta = beta, D = D, n = 500,
#'                 family = family, zeta = c(0, -0.2), sigma = sigma,
#'                 disp.formulas = disp.formulas, gamma = gamma)$data
#'
#' # Now fit using joint
#' long.formulas <- list(
#'   Y.1 ~ time + cont + bin + (1+time|id),
#'   Y.2 ~ time + cont + bin + (1+time|id)
#' )
#' fit <- joint(
#'   long.formulas, Surv(survtime, status) ~ bin,
#'   data, family, disp.formulas = disp.formulas
#' )
#' fit
#' summary(fit)
#' }
joint <- function(long.formulas, surv.formula, 
                  data, family,
                  disp.formulas = NULL, 
                  control = list()){
  
  con <- list(correlated = T, gh.nodes = 3, gh.sigma = 1, center.ph = T,
              tol.abs = 1e-3, tol.rel = 1e-2, tol.thr = 1e-1, tol.den = 1e-3,
              maxit = 200, conv = 'sas', verbose = F, return.inits = F, return.dmats = T, post.process = T)
  conname <- names(con)
  if(any(!names(control)%in%conname)){
    warning("Supplied control arguments do not match with possible names:\n", 
            paste0(sapply(conname, sQuote), collapse=', '), '.')
  }
  con[(conname <- names(control))] <- control
  
  # Ensure supplied families are character, not functions
  if(!is.list(family) | !all(sapply(family, function(x) is.character(x) & !is.function(x))))
    stop(sQuote("family"), " must be supplied as a list of character strings")
  
  start.time <- proc.time()[3]
  
  # Initial parsing ----
  pf <- parent.frame()
  if("factor"%in%class(data$id)) data$id <- as.numeric(as.character(data$id))
  formulas <- lapply(long.formulas, parseFormula)
  surv <- parseCoxph(surv.formula, data, con$center.ph)
  n <- surv$n; K <- length(family)
  if(K!=length(long.formulas)) stop('Mismatched lengths of ', sQuote("family"), " and ", sQuote("long.formulas"),'.')
  # Parse dispersion formulas
  if(is.null(disp.formulas)){
    disp.formulas <- replicate(K, ~1, simplify = FALSE)
  }else{
    if(sum(!sapply(disp.formulas, is.null)) != K)
      stop("Need to supply dispersion formulas for all responses even if not required for all K responses\n",
           "(You can just fill-in ", sQuote("~1"), " for those with no wanted dispersion model.)")
  }
  # Reduce memory overheads, particularly for returned object(s) ?
  disp.formulas <- lapply(disp.formulas, function(x){
    environment(x) <- pf
    x
  })
  
  # Initial conditons ----
  inits.long <- Longit.inits(long.formulas, disp.formulas, data, family)
  id.assign <- AssignIds(data)
  dmats <- getdmats(inits.long$fits, id.assign)
  # Suss out indices of b_k and beta_k.
  b.inds <- lapply(seq_along(dmats$q), function(x){
    xx <- seq(dmats$q[x])
    if(x > 1) 
      return(xx + sum(dmats$q[1:(x-1)]))
    else
      return(xx)
  })
  beta.inds <- lapply(seq_along(dmats$P), function(x){
    xx <- seq(dmats$P[x])
    if(x > 1) 
      return(xx + sum(dmats$P[1:(x-1)]))
    else
      return(xx)
  })
  inds <- list(
    R = list(b = b.inds, beta = beta.inds),
    Cpp = list(b = lapply(b.inds, function(x) x - 1), beta = lapply(beta.inds, function(x) x - 1))
  )
  
  inits.surv <- TimeVarCox(data, inits.long$b, surv, formulas, b.inds, inits.long)
  
  # Longitudinal parameters
  beta <- inits.long$beta.init
  names.beta <- names(beta)
  D <- inits.long$D.init
  sigma <- inits.long$sigma.init # dispersions
  b <- lapply(1:n, function(i) inits.long$b[i, ])
  # Survival parameters
  zeta <- inits.surv$inits[match(names(surv$ph$assign), names(inits.surv$inits))]
  names(zeta) <- paste0('zeta_', names(zeta))
  gamma <- inits.surv$inits[grepl('gamma\\_', names(inits.surv$inits))]
  
  # Survival data objects 
  sv <- surv.mod(surv, formulas, inits.surv$l0.init, inits.long)
  
  # Parameter vector and list ----
  Omega <- list(D=D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
  params <- c(setNames(vech(D), paste0('D[', apply(which(lower.tri(D, T), arr.ind = T), 1, paste, collapse = ','), ']')),
              beta, unlist(sigma)[inits.long$sigma.include], gamma, zeta)
  sigma.include <- inits.long$sigma.include
  if(!con$return.inits) rm(inits.surv)
  
  # Gauss-Hermite Quadrature ----
  GH <- statmod::gauss.quad.prob(con$gh.nodes, 'normal', sigma = con$gh.sigma)
  w <- GH$w; v <- GH$n

  # Begin EM ----
  diff <- 100; iter <- 0;
  # Convergence criteria setup
  if(!con$conv%in%c('absolute', 'relative', 'either', 'sas')){
    warning('Convergence criteria must be one of "absolute", "relative", "either" or "sas". Using "sas"')
    con$conv <- "sas"
  }
    
  convergence.criteria <- list(type = con$conv, tol.abs = con$tol.abs, tol.rel = con$tol.rel, tol.den = con$tol.den,
                               threshold = con$tol.thr)
  
  if(con$verbose){
    cat("Initial conditions: \n")
    converge.check(0, 0, convergence.criteria, 0, Omega, TRUE)
    cat("\nStarting EM algorithm...\n")
  }
  converged <- FALSE
  EMstart <- proc.time()[3]
  while((!converged) && (iter < con$maxit)){
    update <- EMupdate(Omega, family, dmats, b, sv, 
                       surv, w, v, con, inds)
    
    if(!con$correlated) update$D[inits.long$off.inds] <- 0
    params.new <- c(vech(update$D), update$beta, unlist(update$sigma)[sigma.include], 
                    update$gamma, update$zeta)
    names(params.new) <- names(params)
    
    # Set new estimates as current
    b <- update$b
    D <- update$D; beta <- update$beta; sigma <- update$sigma
    gamma <- update$gamma; zeta <- update$zeta
    l0 <- update$l0; l0u <- update$l0u; l0i <- update$l0i
    sv$l0 <- l0; sv$l0i <- l0i; sv$l0u <- l0u
    iter <- iter + 1
    Omega <- list(D = D, beta = beta, sigma = sigma, gamma = gamma, zeta = zeta)
    
    convcheck <- converge.check(params, params.new, convergence.criteria, iter, Omega, con$verbose)
    if(iter >= 4) converged <- convcheck$converged # Allow to converge after 3 iterations
    params <- params.new
  }
  
  EMend <- proc.time()[3]
  coeffs <- Omega
  coeffs$beta <- setNames(c(Omega$beta), names.beta)
  if(is.not.SPD(coeffs$D)) warning("Covariance matrix D is not positive semi-definite at convergence,",
                                   " potential model misspecification or lower tolerance options required.")
  
  out <- list(coeffs = coeffs,
              hazard = cbind(ft = sv$ft, haz = l0, nev = sv$nev))
  
  ModelInfo <- list()
  ModelInfo$ResponseInfo <- sapply(1:K, function(k){
    paste0(inits.long$responses[k], ' (', family[k], ')')
  })
  ModelInfo$Resps <- inits.long$responses
  ModelInfo$family <- family
  ModelInfo$K <- dmats$K
  ModelInfo$Pcounts <- list(P = dmats$P, Pd = dmats$Pd, 
                            q = sv$q, vD = length(vech(D)))
  ModelInfo$long.formulas <- long.formulas
  ModelInfo$disp.formulas <- disp.formulas
  ModelInfo$surv.formula <- surv.formula
  ModelInfo$survtime <- surv$survtime
  ModelInfo$status <- surv$status
  ModelInfo$control <- con
  ModelInfo$convergence.criteria <- convergence.criteria
  ModelInfo$inds <- inds
  ModelInfo$n <- dmats$n
  ModelInfo$nobs <- setNames(dmats$m, inits.long$responses)
  ModelInfo$mi <- sapply(dmats$mi, unlist)
  ModelInfo$nev <- sum(sv$nev)
  ModelInfo$id.assign <- list(original.id = id.assign$id,
                              assigned.id = id.assign$assign)
  out$ModelInfo <- ModelInfo
  
  # Post processing ----
  if(con$post.process){
    if(con$verbose) cat('Post-processing...\n')
    gamma.rep <- rep(gamma, sapply(b.inds, length))
    pp.start.time <- proc.time()[3]
    
    II <- obs.emp.I(coeffs, dmats, surv, sv, family, b, 
                    l0i, l0u, w, v, inds)
    H <- structure(II$Hessian,
                   dimnames = list(names(params), names(params)))
    
    I <- tryCatch(solve(H), error = function(e) e)
    if(inherits(I, 'error')) I <- structure(MASS::ginv(H), dimnames = dimnames(H))
    out$Hessian <- H
    out$vcov <- I
    out$SE <- sqrt(diag(I))
   
    postprocess.time <- round(proc.time()[3] - pp.start.time, 2)
    
    # Calculate log-likelihood. Done separately as EMtime + postprocess.time is for EM + SEs.
    out$logLik <- joint.log.lik(coeffs, dmats, surv, sv, family, II$b.hat, l0i, l0u, inds, II$Sigma)
    # Collate RE and their variance
    REs <- do.call(rbind, II$b.hat)
    attr(REs, 'Var') <- do.call(rbind, lapply(II$Sigma, diag))
    attr(REs, 'vcov') <- do.call(rbind, lapply(II$Sigma, vech))
    out$REs <- REs
  }else{
    REs <- do.call(rbind, b)
    attr(REs, 'Var') <- do.call(rbind, lapply(update$Sigma, diag))
    attr(REs, 'vcov') <- do.call(rbind, lapply(update$Sigma, vech))
    out$REs <- REs
  }
  comp.time <- round(proc.time()[3] - start.time, 3)
  out$elapsed.time <- c(`EM time` = unname(round(EMend - EMstart, 3)),
                        `Post processing` = if(con$post.process) unname(postprocess.time) else NULL,
                        `Total Computation time` = unname(comp.time),
                        `iterations` = iter)
  if(con$post.process) sv <- surv.mod(surv, formulas, l0, inits.long)
  dmats <- list(long = dmats, surv = sv, ph = surv)
  if(con$return.dmats) out$dmats <- dmats
  
  if(con$return.inits) out$inits = list(inits.long = inits.long,
                                        inits.surv = inits.surv)
  class(out) <- 'joint'
  return(out)
}

#' @method print joint
#' @keywords internal
#' @export
print.joint <- function(x, ...){
  if(!inherits(x, 'joint')) stop('x must be a "joint" object!')
  
  M <- x$ModelInfo
  K <- M$K # Number of responses
  fams <- M$family
  dpsL <- sapply(M$long.formulas, long.formula.to.print, 1)
  dpsD <- sapply(M$disp.formulas, deparse)
  dpsS <- deparse(M$surv.formula)
  
  # Data information
  cat(sprintf("Number of subjects: %d\n", M$n))
  cat(sprintf("Number of events: %d (%.2f%%)\n", M$nev, 100 * M$nev/M$n))
  
  # Longitudinal information: 
  cat("\n===================\nModel specification\n===================\n")
  if(K == 1)
    cat("Univariate longitudinal process specification:\n")
  else
    cat("Multivariate longitudinal process specifications: \n")
  
  for(k in 1:K){
    cat(sprintf("%s: %s\n", M$ResponseInfo[k], dpsL[k]))
    if(fams[[k]] %in% c("Gamma", "negbin", "genpois"))
      cat(sprintf("Dispersion model: %s\n",
                  if(dpsD[k]=="~1") "(Intercept) model" else dpsD[k]
      ))
  }
  
  cat("\nSurvival sub-model specification: \n")
  cat(deparse(M$surv.formula))
  
  cat("\n\nAssociation parameter estimates: \n")
  print(setNames(x$coeffs$gamma,
        M$Resps))
  cat("\n")
  invisible(x)
}
