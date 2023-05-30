#' Simulate realisations from a generalised poisson distribution
#' @param mu A numeric vector of rates \eqn{\exp{\eta}}, with \eqn{\eta} the linear predictor.
#' @param phi A numeric specifying the dispersion \eqn{\varphi}. If \eqn{\varphi<0} the response 
#' will be under-dispersed and overdispersed if \eqn{\varphi>0}.
#' 
#' @details Follows the "GP-1" implementation of the generalised Poisson distribution outlined 
#' in Zamani & Ismail (2012). The variance of produced \eqn{Y} is \eqn{(1+\varphi)^2\mu}.
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
    kum <- GP1_pmf_scalar(mu[j], phi, 0)
    while(rand > kum){
      ans <- ans + 1
      # message(ans)
      kum <- kum + GP1_pmf_scalar(mu[j], phi, ans)
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
#' @param family a \eqn{K}-list of families, see \strong{details}.
#' @param sigma a \eqn{K}-vector of dispersion parameters corresponding to the order of 
#' \code{family}; see \strong{details}.
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
#' @param unif.times logical, if \code{unif.times = TRUE} (the default), then \emph{every} subject
#' will have the same follow-up times defined by \code{seq(0, fup, length.out = ntms)}. If 
#' \code{unif.times = FALSE} then follow-up times are set as random draws from a uniform 
#' distribution with maximum \code{fup}. 
#' @param dof integer, specifies the degrees of freedom of the multivariate t-distribution
#' used to generate the random effects. If specified, this t-distribution is used. If left
#' at the default \code{dof=Inf} then the random effects are drawn from a multivariate normal
#' distribution.
#' @param random.formula allows user to specify if an intercept-and-slope (\code{~ time}) or 
#' intercept-only (\code{~1}) random effects structure should be used. defaults to the former.
#' @param return.ranefs a logical determining whether the \emph{true} random effects should be 
#' returned. This is largely for internal/simulation use. Default \code{return.ranefs = FALSE}.
#'   
#' @returns A list of two \code{data.frame}s: One with the full longitudinal data, and another 
#' with only survival data. If \code{return.ranefs=TRUE}, a matrix of the true \eqn{b} values is
#' also returned.
#'
#' @details \code{simData} simulates data from a multivariate joint model with a mixture of 
#' families for each \eqn{K=1,\dots,3} response. Currently, the argument \code{random.formula}
#' specifies the association structure for \strong{all} responses. The specification of 
#' \code{family} changes requisite dispersion parameter, if applicable. The \code{family} list can
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
#'   \item{\code{"genpois"}}{Simulated with a log link, corresponding item in \code{sigma} will be
#'   the \strong{dispersion}. Values < 0 correspond to under-dispersion, and values > 0 over-
#'   dispersion. See \code{\link{rgenpois}} for more information. Simulated variance is 
#'   \eqn{(1+\varphi)^2\mu}.}
#'   \item{\code{"Gamma"}}{Simulated with a log link, corresponding item in \code{sigma} will be
#'   the \strong{shape}.}
#'   
#'   }
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
#' 
#' # K = 3 mixture of families with dispersion parameters
#' beta <- do.call(rbind, replicate(3, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
#' gamma <- c(0.3, -0.3, 0.3)
#' D <- diag(c(0.25, 0.09, 0.25, 0.05, 0.25, 0.09))
#' family <- list('gaussian', 'genpois', 'Gamma')
#' sigma <- c(.16, 1.5, 1.5)
#' sim.data <- simData(ntms=15, family = family, sigma = sigma, beta = beta, D = D, gamma = gamma,
#'                     theta = c(-3, 0.2), zeta = c(0,-.2))
#' 
#' # K = 4 mixture of families with/out dispersion
#' beta <- do.call(rbind, replicate(4, c(2, -0.1, 0.1, -0.2), simplify = FALSE))
#' gamma <- c(-0.75, 0.3, -0.6, 0.5)
#' D <- diag(c(0.25, 0.09, 0.25, 0.05, 0.25, 0.09, 0.16, 0.02))
#' family <- list('gaussian', 'poisson', 'binomial', 'gaussian')
#' sigma <- c(.16, 0, 0, .05)
#' sim.data <- simData(ntms=15, family = family, sigma = sigma, beta = beta, D = D, gamma = gamma,
#'                     theta = c(-3, 0.2), zeta = c(0,-.2))
simData <- function(n = 250, ntms = 10, fup = 5, 
                    family = list('gaussian', 'gaussian'), 
                    sigma = c(0.16, 0.16),
                    beta = rbind(c(1, 0.10, 0.33, -0.50), c(1, 0.10, 0.33, -0.50)), D = NULL,
                    gamma = c(0.5, -0.5), zeta = c(0.05, -0.30),
                    theta = c(-4, 0.2), cens.rate = exp(-3.5),
                    unif.times = TRUE,
                    dof = Inf,
                    random.formula = NULL,
                    return.ranefs = FALSE){
  
  # Checks --------------
  # Check family is valid option
  family <- lapply(family, function(f){
    if("function"%in%class(f)) f <- f()$family
    if(!f%in%c("gaussian", "binomial", "poisson", "genpois", "Gamma")) 
      stop('Family must be one of "gaussian", "binomial", "poisson", "genpois" or "Gamma".')
    f
  })
  # Check dispersion supplied if neg. binom.
  family <- lapply(family, function(f) family <- match.arg(f, c("gaussian", "binomial", "poisson", "genpois", "Gamma"), several.ok = F))
  
  # Checks wrt the fixed effects and dispersion parameters
  funlist <- unlist(family)
  num.gp <- length(which(funlist == 'genpois'))
  num.ga <- length(which(funlist == 'gaussian'))
  num.Ga <- length(which(funlist == 'Gamma'))
  
  # Fixed effects
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
  if(any(eigen(D)$value < 0) || det(D) <= 0) stop('D must be positive semi-definite')
  
  # Necessary parameters & data generation ----
  if(unif.times){
    time <- seq(0, fup, length.out = ntms)
    tau <- fup + 0.1
  }else{
    time <- c(replicate(n, c(0, sort(runif(ntms - 1, max = fup)))))
    tau <- fup
  }

  cont <- rnorm(n); bin <- rbinom(n, 1, 0.5)  # Continuous and binary covariates.
  
  df <- data.frame(id = rep(1:n, each = ntms),
                   time = if(unif.times) rep(time, n) else time,
                   cont = rep(cont, each = ntms),
                   bin = rep(bin, each = ntms))
  
  # Design matrices, used across ALL responses ----
  X <- model.matrix(~ time + cont + bin, data = df)
  if(!is.null(random.formula)){
    Z <- lapply(random.formula, function(x) model.matrix(x, data = df))
  }else{
    Z <- replicate(K, model.matrix(~ 1 + time, data = df), simplify = F) # assume all intslopes.
  }

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
    switch(f, 
           gaussian = Y <- rnorm(n * ntms, mean = etak, sd = sqrt(sigma[k])),
           binomial = Y <- rbinom(n * ntms, 1, plogis(etak)),
           poisson = Y <- rpois(n * ntms, exp(etak)),
           Gamma = Y <- rgamma(n * ntms, shape = sigma[k], scale = exp(etak)/sigma[k]),
           genpois = {
             maxtry <- 5
             try <- 1
             e <- T
             while(e){ # genpois can be quite sensitive in simulation part; repeat sim up to 5 times.
               Y <- tryCatch(rgenpois(exp(etak), sigma[k]),
                             error = function(e) NA)
               e <- any(is.na(Y))
               try <- try + 1
               if(try > maxtry) stop("Issues creating genpois response with beta values ", c(betak), 
                                     " and mean(b)", c(round(mean(b[df$id, b.inds[[k]]]), 3)))
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
  message(round(100 * sum(surv.data$status)/n), '% failure rate')
  
  out <- list(data =  out.data, 
              surv.data =  out.data[!duplicated(out.data[,'id']), c('id', 'survtime', 'status', 'cont', 'bin')])
  if(return.ranefs) out$ranefs <- b
  out
}
