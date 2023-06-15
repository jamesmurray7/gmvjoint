#' @keywords internal
dgenpois <- function(x, mu, phi, log = FALSE){
  lx <- length(x); lm <- length(mu); lp <- length(phi)
  ch <- lx == lm & lx == lp & lm == lp
  if(!ch) stop("Not all supplied inputs are of same length, check.\n")
  
  out <- sapply(1:lx, function(i) GP1_pmf_scalar(mu[i], phi[i], x[i]))
  if(log)
    return(log(out))
  else
    return(out)
}

#' Simulate realisations from a generalised poisson distribution
#' @param mu A numeric vector of rates \eqn{\exp{\eta}}, with \eqn{\eta} the linear predictor.
#' @param phi A numeric specifying the dispersion \eqn{\varphi}. If \eqn{\varphi<0} the response 
#' will be under-dispersed and overdispersed if \eqn{\varphi>0}.
#' 
#' @details Follows the "GP-1" implementation of the generalised Poisson distribution outlined 
#' in Zamani & Ismail (2012). The variance of produced \eqn{Y} is \eqn{(1+\varphi)^2\mu}. As such
#' the dispersion parameter is bounded (i.e. not in positive reals as with CMP distribution).
#' 
#' @return An appropriately-dimensioned vector of count data.
#' 
#' @references 
#' 
#' Zamani H and Ismail N. Functional Form for the Generalized Poisson Regression Model, 
#' \emph{Communications in Statistics - Theory and Methods} 2012; \strong{41(20)}; 3666-3675.
#'
#' @export
rgenpois <- function(mu, phi){
  n <- length(mu)
  out <- numeric(n)
  for(j in 1:n){
    ans <- 0;
    rand <- runif(1)
    kum <- dgenpois(0, mu[j], phi[j])
    while(rand > kum){
      ans <- ans + 1
      kum <- kum + dgenpois(ans, mu[j], phi[j])
      # For debugging
      if(is.infinite(kum) || is.nan(kum)) cat(sprintf("NaN/Infinite kum: mu = %.2f, phi = %.2f", mu[j], phi[j]))
    }
    out[j] <- ans
  }
  out
}

#' Simulate data from a multivariate joint model
#' 
#' @description Simulate multivariate longitudinal and survival data from a joint model 
#' specification, with potential mixture of response families. Implementation is similar 
#' to existing packages (e.g. \code{joineR}, \code{joineRML}).
#'   
#' @param n the number of subjects
#' @param ntms the number of time points
#' @param fup the maximum follow-up time, such that t = [0, ..., fup] with length \code{ntms}. 
#' In instances where subject \eqn{i} \emph{doesn't} fail before \code{fup}, their censoring
#' time is set as \code{fup + 0.1}.
#' @param family a \eqn{K}-\code{list} of families, see \strong{details}.
#' @param sigma a \eqn{K}-\code{list} of dispersion parameters corresponding to the order of 
#' \code{family}, and matching \code{disp.formulas} specification; see \strong{details}.
#' @param beta a \eqn{K \times 4} matrix specifying fixed effects for each \eqn{K} parameter, 
#' in the order (Intercept), time, continuous, binary.
#' @param D a positive-definite matrix specifying the variance-covariance matrix for the random
#' effects. If not supplied an identity matrix is assumed.
#' @param gamma a \eqn{K}-vector specifying the association parameters for each longitudinal 
#' outcome.
#' @param zeta a vector of length 2 specifying the coefficients for the baseline covariates in 
#' the survival sub-model, in the order of continuous and binary.
#' @param theta parameters to control the failure rate, see \strong{baseline hazard}.
#' @param cens.rate parameter for \code{rexp} to generate censoring times for each subject.
#' @param regular.times logical, if \code{regular.times = TRUE} (the default), then 
#' \emph{every} subject will have the same follow-up times defined by 
#' \code{seq(0, fup, length.out = ntms)}. If \code{regular.times = FALSE} then follow-up times are
#' set as random draws from a uniform distribution with maximum \code{fup}. 
#' @param dof integer, specifies the degrees of freedom of the multivariate t-distribution
#' used to generate the random effects. If specified, this t-distribution is used. If left
#' at the default \code{dof=Inf} then the random effects are drawn from a multivariate normal
#' distribution.
#' @param random.formulas allows user to specify if an intercept-and-slope (\code{~ time}) or 
#' intercept-only (\code{~1}) random effects structure should be used on a response-by-response
#' basis. Defaults to an intercept-and-slope for all responses.
#' @param disp.formulas allows user to specify the dispersion model simulated. Intended use is
#' to allow swapping between intercept only (the default) and a time-varying one (\code{~ time}).
#' Note that this should be a \eqn{K}-\code{list} of formula objects, so if only one dispersion 
#' model is wanted, then an intercept-only should be specified for remaining sub-models. The
#' corresponding item in list of \code{sigma} parameters should be of appropriate size. Defaults
#' to an intercept-only model.
#' @param return.ranefs a logical determining whether the \emph{true} random effects should be 
#' returned. This is largely for internal/simulation use. Default \code{return.ranefs = FALSE}.
#'
#' @returns A list of two \code{data.frame}s: One with the full longitudinal data, and another 
#' with only survival data. If \code{return.ranefs=TRUE}, a matrix of the true \eqn{b} values is
#' also returned. By default (i.e. no arguments provided), a bivariate Gaussian set of joint
#' data is returned.
#'
#' @details \code{simData} simulates data from a multivariate joint model with a mixture of 
#' families for each \eqn{k=1,\dots,K} response. The specification of \code{family} changes
#' requisite dispersion parameter \code{sigma}, if applicable. The \code{family} list can
#' (currently) contain: 
#'   
#'   \describe{
#'  
#'   \item{\code{"gaussian"}}{Simulated with identity link, corresponding item in \code{sigma}
#'   will be the \strong{variance}.}
#'   \item{\code{"poisson"}}{Simulated with log link, corresponding dispersion in \code{sigma} 
#'   can be anything, as it doesn't impact simulation.}
#'   \item{\code{"binomial"}}{Simulated with logit link, corresponding dispersion in \code{sigma} 
#'   can be anything, as it doesn't impact simulation.}
#'   \item{\code{"negbin"}}{Simulated with a log link, corresponding item in \code{sigma} will be
#'   the \strong{overdispersion} defined on the log scale. Simulated variance is 
#'   \eqn{\mu+\mu^2/\varphi}.}
#'   \item{\code{"genpois"}}{Simulated with a log link, corresponding item in \code{sigma} will be
#'   the \strong{dispersion}. Values < 0 correspond to under-dispersion, and values > 0 over-
#'   dispersion. See \code{\link{rgenpois}} for more information. Simulated variance is 
#'   \eqn{(1+\varphi)^2\mu}.}
#'   \item{\code{"Gamma"}}{Simulated with a log link, corresponding item in \code{sigma} will be
#'   the \strong{shape} parameter, defined on the log-scale.}
#'   
#'   }
#'   
#'  Therefore, for families \code{"negbin"}, \code{"Gamma"}, \code{"genpois"}, the user can
#'  define the dispersion model desired in \code{disp.formulas}, which creates a data matrix 
#'  \eqn{W}. For the \code{"negbin"} and \code{"Gamma"} cases, we define 
#'  \eqn{\varphi_i=\exp\{W_i\sigma_i\}} (i.e. the exponent of the linear predictor of the 
#'  dispersion model); and for \code{"genpois"} the identity of the linear is used.
#'   
#' @section Baseline hazard: 
#'  
#'   When simulating the survival time, the baseline hazard is a Gompertz distribution controlled 
#'   by \code{theta=c(x,y)}:
#'   
#'   \deqn{\lambda_0(t) = \exp{x + yt}}
#'   
#'   where \eqn{y} is the shape parameter, and the scale parameter is \eqn{\exp{x}}.
#'
#' @importFrom MASS mvrnorm
#' @importFrom mvtnorm rmvt
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @keywords simulation
#' @seealso \code{\link{joint}}
#' 
#' @references 
#' 
#' Austin PC. Generating survival times to simulate Cox proportional hazards
#' models with time-varying covariates. \emph{Stat Med.} 2012; \strong{31(29)}:
#' 3946-3958.
#' 
#' @export 
#'
#' @examples
#' # 1) A set of univariate data ------------------------------------------
#' beta <- c(2.0, 0.33, -0.25, 0.15)
#' # Note that by default arguments are bivariate, so need to specify the univariate case
#' univ.data <- simData(beta = beta,    
#'                      gamma = 0.15, sigma = list(0.2), family = list("gaussian"), 
#'                      D = diag(c(0.25, 0.05)))
#'                      
#' # 2) Univariate data, with failure rate controlled ---------------------
#' # In reality, many facets contribute to the simulated failure rate, in 
#' # this example, we'll just atler the baseline hazard via 'theta'.
#' univ.data.highfail <- simData(beta = beta,
#'                               gamma = 0.15, sigma = list(0.0), family = list("poisson"),
#'                               D = diag(c(0.40, 0.08)), theta = c(-2, 0.1))
#' 
#' # 3) Trivariate (K = 3) mixture of families with dispersion parameters -
#' beta <- do.call(rbind, replicate(3, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
#' gamma <- c(0.3, -0.3, 0.3)
#' D <- diag(c(0.25, 0.09, 0.25, 0.05, 0.25, 0.09))
#' family <- list('gaussian', 'genpois', 'negbin')
#' sigma <- list(.16, 1.5, log(1.5))
#' triv.data <- simData(ntms=15, family = family, sigma = sigma, beta = beta, D = D, 
#'                      gamma = gamma, theta = c(-3, 0.2), zeta = c(0,-.2))
#' 
#' # 4) K = 4 mixture of families with/out dispersion ---------------------
#' beta <- do.call(rbind, replicate(4, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
#' gamma <- c(-0.75, 0.3, -0.6, 0.5)
#' D <- diag(c(0.25, 0.09, 0.25, 0.05, 0.25, 0.09, 0.16, 0.02))
#' family <- list('gaussian', 'poisson', 'binomial', 'gaussian')
#' sigma <- list(.16, 0, 0, .05) # 0 can be anything here, as it is ignored internally.
#' mix.data <- simData(ntms=15, family = family, sigma = sigma, beta = beta, D = D, gamma = gamma,
#'                     theta = c(-3, 0.2), zeta = c(0,-.2))
#'                     
#' # 5) Bivariate joint model with two dispersion models. -----------------
#' disp.formulas <- list(~time, ~time)          # Two time-varying dispersion models
#' sigma <- list(c(0.00, -0.10), c(0.10, 0.15)) # specified in form of intercept, slope
#' D <- diag(c(.25, 0.04, 0.50, 0.10))
#' disp.data <- simData(family = list("genpois", "negbin"), sigma = sigma, D = D,
#'                      beta = rbind(c(0, 0.05, -0.15, 0.00), 1 + c(0, 0.25, 0.15, -0.20)),
#'                      gamma = c(1.5, 1.5),
#'                      disp.formulas = disp.formulas, fup = 5)            
#'
#' # 6) Trivariate joint model with mixture of random effects models ------
#' # It can be hard to e.g. fit a binomial model on an intercept and slope, since e.g.
#' # glmmTMB might struggle to accurately fit it (singular fits, etc.). To that end, could
#' # swap the corresponding random effects specification to be an intercept-only.
#' family <- list("gaussian", "binomial", "gaussian")
#' # A list of formulae, even though we want to change the second sub-model's specification
#' # we need to specify the rest of the items, too (same as disp.formulas, sigma).
#' random.formulas <- list(~time, ~1, ~time)
#' beta <- rbind(c(2, -0.2, 0.5, -0.25), c(0, 0.5, 1, -1), c(-2, 0.2, -0.5, 0.25))
#' # NOTE that the specification of RE matrix will need to match.
#' D <- diag(c(0.25, 0.09, 1, 0.33, 0.05))
#' # Simulate data, and return REs as a sanity check...
#' mix.REspec.data <- simData(beta = beta, D = D, family = family,
#'                            gamma = c(-0.5, 1, 0.5), sigma = list(0.15, 0, 0.15),
#'                            random.formulas = random.formulas, return.ranefs = TRUE)
simData <- function(n = 250, ntms = 10, fup = 5, 
                    family = list('gaussian', 'gaussian'), 
                    sigma = list(0.16, 0.16),
                    beta = rbind(c(1, 0.10, 0.33, -0.50), c(1, 0.10, 0.33, -0.50)), D = NULL,
                    gamma = c(0.5, -0.5), zeta = c(0.05, -0.30),
                    theta = c(-4, 0.2), cens.rate = exp(-3.5),
                    regular.times = TRUE,
                    dof = Inf,
                    random.formulas = NULL, disp.formulas = NULL,
                    return.ranefs = FALSE){
  random.formula <- random.formulas
  # Checks --------------
  # Check family is valid option
  family <- lapply(family, function(f){
    if("function"%in%class(f)) f <- f()$family
    if(!f%in%c('gaussian', 'binomial', 'poisson', 'genpois', 'Gamma', 'genpois', 'negbin')) 
      stop('Family must be one of "gaussian", "binomial", "poisson", "genpois", "negbin" or "Gamma".')
    f
  })
  # Check dispersion supplied if neg. binom.
  family <- lapply(family, function(f) match.arg(f, c('gaussian', 'binomial', 'poisson', 'genpois', 'Gamma', 'genpois', 'negbin'), several.ok = F))
  funlist <- unlist(family)
  # Fixed effects
  if(!"matrix"%in%class(beta)) beta <- t(beta) # Assume if a vector is passed as beta, then a univariate set of data is being attempted.
  K <- nrow(beta)
  if(K != length(funlist)) stop('Incorrect number of families provided and/or beta terms.')
  if(K != length(gamma)) stop('Incorrect number of association parameters provided wrt beta terms.')
  if(is.null(D)) D <- diag(K * 2) else D <- D
  if(!is.null(random.formula)){
    if(length(random.formula) != K) stop('Provided random formulas must be a list of length K.')
  }else{
    if((K*2) != ncol(D)) stop('Incorrect dimension on supplied covariance matrix D, do you need to supply random.formula?')
  }
  
  # Check covariance matrix D is positive semi-definite.
  if(is.not.SPD(D)) stop('D must be positive semi-definite')
  
  # Parse dispersion formulas
  pf <- parent.frame()
  if(is.null(disp.formulas)){
    disp.formulas <- replicate(K, ~1, simplify = FALSE)
    lapply(disp.formulas, function(x) environment(x) <- pf)
  }else{
    if(sum(!sapply(disp.formulas, is.null)) != K)
      stop("Need to supply dispersion formulas for all responses even if not required for all K responses\n",
           "(You can just fill-in ~1 for those with no wanted dispersion model.)")
  }
  
  # Necessary parameters & data generation ----
  if(regular.times){
    time <- seq(0, fup, length.out = ntms)
    tau <- fup + 0.1
  }else{
    time <- c(replicate(n, c(0, sort(runif(ntms - 1, max = fup)))))
    tau <- fup
  }

  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = if(regular.times) rep(time, n) else time,
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices, used across ALL responses ----
  X <- model.matrix(~ time + cont + bin, data = df)
  W <- lapply(disp.formulas, function(x) model.matrix(x, data = df))
  if(!is.null(random.formula)){
    Z <- lapply(random.formula, function(x) model.matrix(x, data = df))
  }else{
    Z <- replicate(K, model.matrix(~ 1 + time, data = df), simplify = F) # assume all intslopes.
  }
  # Check that the supplied disp.formulas match with elements of sigma
  if(!all(sapply(1:K, function(k) ncol(W[[k]])==length(sigma[[k]]))))
    stop("Mis-specified sigma or disp.formulas.\n")

  # Random effects specification
  if(is.infinite(dof)) 
    .simRE <- function(n, mu, Sigma) MASS::mvrnorm(n = n, mu = mu, Sigma = Sigma)
  else
    .simRE <- function(n, mu, Sigma) mvtnorm::rmvt(n = n, sigma = Sigma, delta = mu, df = dof)
  
  # Linear predictor & response generation ----
  if(!is.null(random.formula)){
    q <- ncol(do.call(cbind, lapply(Z, head)))
    # b <- MASS::mvrnorm(n, mu = rep(0, q), Sigma = D)
    b <- .simRE(n, rep(0, q), D)
    b.inds <- split(1:q, do.call(c, sapply(1:K, function(k) rep(k, ncol(Z[[k]])), simplify = F)))
  }else{
    # b <- MASS::mvrnorm(n, mu = rep(0, K * 2), Sigma = D)
    b <- .simRE(n, rep(0, K * 2), D)
    b.inds <- split(1:(2*K), rep(1:K, each = 2)) 
  }
  
  Y <- sapply(1:K, function(k){
    f <- family[[k]]; Zk <- Z[[k]]
    betak <- beta[k,,drop=T]
    etak <- X %*% betak + rowSums(Zk * b[df$id, b.inds[[k]]])
    phik <- W[[k]] %*% sigma[[k]]
    switch(f, 
           gaussian = Y <- rnorm(n * ntms, mean = etak, sd = sqrt(phik[1])),
           binomial = Y <- rbinom(n * ntms, 1, plogis(etak)),
           poisson = Y <- rpois(n * ntms, exp(etak)),
           Gamma = Y <- rgamma(n * ntms, shape = exp(phik), scale = exp(etak)/exp(phik)),
           negbin = Y <- MASS::rnegbin(n * ntms, mu = exp(etak), theta = exp(phik)),
           genpois = {
             maxtry <- 5
             try <- 1
             e <- T
             while(e){ # genpois can be quite sensitive in simulation part; repeat sim up to 5 times.
               Y <- tryCatch(rgenpois(exp(etak), phik),
                             error = function(e) NA)
               e <- any(is.na(Y))
               try <- try + 1
               if(try > maxtry){
                 stop(sprintf("Issues creating genpois response, this could be due to large values for mu (range: %.2f-%.2f),",
                              range(exp(etak))[1], range(exp(etak))[2]),
                      sprintf(" or dispersion (range: %.2f, %.2f); it's recommended that this value is between", 
                              range(phik)[1], range(phik)[2]),
                      sprintf(" %.2f and %.2f.\n", max(-1, -max(exp(etak)/4)), 1))
               }
             }
             Y
           })
  })
  colnames(Y) <- paste0('Y.', 1:K)
  
  df <- cbind(df, Y)
  
  # Survival ----
  theta0 <- theta[1]; theta1 <- theta[2]
  Keta <- cbind(cont, bin) %*% zeta
  U <- runif(n)
  if(!is.null(random.formula)){
    zz <- lapply(Z, colnames)
    ints <- which(grepl('\\(Intercept\\)', do.call(c, zz))); slopes <- which(grepl('time', do.call(c, zz)))
    if(any(!grepl('\\(Intercept\\)|time', do.call(c, zz)))) message('Warning: Only intercept-only or intercept-ands-slope random.formulas allowed.')
    # Which of the K responses have intercept/slope
    ints.gamma <- do.call(c, lapply(1:K, function(k){
      if(any(grepl('\\(Intercept\\)', colnames(Z[[k]])))) return(k) else return(NULL)
    }))
    slopes.gamma <- do.call(c, lapply(1:K, function(k){
      if(any(grepl('time', colnames(Z[[k]])))) return(k) else return(NULL)
    }))
  }else{
    ints <- seq(1,2*K,by=2); slopes <- seq(2, 2*K, by = 2) # intslope on all.
    ints.gamma <- slopes.gamma <- 1:K    # Each K has both intercept and slope.
  }
  
  b0 <- b[, ints, drop = F]; b1 <- b[, slopes, drop = F]
  # Generate survival times (Austin et al 2012)
  denom <- theta1 + b1 %*% gamma[slopes.gamma]
  rhs <- (theta1 + b1 %*% gamma[slopes.gamma]) * log(U)/(exp(theta0 + Keta + b0 %*% gamma[ints.gamma]))
  t <- suppressWarnings(log(1-rhs)/denom)
  t[is.nan(t)] <- tau
  
  # Collect survival times, and generate cenors times.
  cens.time <- rexp(n, cens.rate)
  survtime <- pmin(t, cens.time)
  survtime[survtime >= tau] <- tau
  
  # Status flag
  status <- rep(1, n)
  is.censored <- cens.time < survtime
  status[which(survtime == tau | is.censored | survtime == cens.time)] <- 0 # Failure flag
  
  # Output Dataframes ----
  surv.data <- data.frame(id = 1:n, survtime, status)
  long.data <- df
  
  out.data <- merge(df, surv.data, by = 'id')
  out.data <- out.data[out.data$time < out.data$survtime, ]
  cat(sprintf("%.2f%% failure rate.\n", 100 * sum(surv.data$status)/n))
  
  out <- list(data =  out.data, 
              surv.data =  out.data[!duplicated(out.data[,'id']), c('id', 'survtime', 'status', 'cont', 'bin')])
  if(return.ranefs) out$ranefs <- b
  out
}
