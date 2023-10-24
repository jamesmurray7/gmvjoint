#' Dynamic predictions for survival sub-model in a multivariate joint model.
#'
#' @description Calculates individualised conditional survival probabilities for subjects
#' during a period of follow-up using a \code{joint} model fit along with requisite longitudinal 
#' process history. 
#' 
#' \strong{Note} that this function is largely designed for use within the ROC function
#' which assesses discriminatory power of the joint model, however it \emph{does} function
#' by itself with proper use of its arguments.
#' 
#' @param data the data to which the original \code{joint} model was fit.
#' @param id subject identifier, i.e. for which subject is the conditional survival probabilities 
#' desired?
#' @param fit a joint model fit by the  \code{joint} function. 
#' @param u a numeric \code{vector} of candidate follow-up times for which a survival probability 
#' should be calculated. Note that the first item \code{u[1]} denotes the start of the "window"
#' and is dropped from calculations. If \code{u=NULL} (the default), then the probability of 
#' surviving all failure times after the \code{id}'s final longitudinal \code{time} is calculated.
#' @param nsim how many Monte Carlo simulations should be carried out? Defaults to 
#' \code{nsim=200}. First-order estimates are calculated if \code{nsim=0}.
#' @param progress a logical, if \code{progress=TRUE} (the default) then a progress bar displaying
#' the current percentage of simulations have been completed.
#' @param scale numeric scales the variance-covariance parameter in the proposal distribution for 
#' the Metropolis-Hastings algorithm. Defaults to \code{scale = NULL} which doesn't scale the
#' variance term at all. Users are encouraged to experiment with values here; this parameter
#' controls the acceptance rate of the MH scheme.
#' @param df numeric denotes the degrees of freedom of the proposed \eqn{t} distribution on
#' the random effects; \code{df=4} is suggested and is the default.
#'
#' @return A list of class \code{dynPred} which consists of three items: \describe{
#'   \item{\code{pi}}{A \code{data.frame} which contains each candidate failure time (supplied by
#'   \code{u}), with the mean, median and 2.5\% and 97.5\% quantiles of probability of survival 
#'   until this failure time.}
#'   \item{\code{pi.raw}}{A \code{matrix} of with \code{nsim} rows and \code{length(u)} columns,
#'   each row represents the \eqn{l}th conditional survival probability of survival each \code{u}
#'   survival time. This is largely for debugging purposes.}
#'   \item{MH.accept}{The acceptance rate of the Metropolis-Hastings algorithm on the random
#'   effects.}
#' }
#' 
#' @details Dynamic predictions for the time-to-event process based on information available
#' on the subject's longitudinal process up to given time \eqn{t} are calculated by Monte Carlo
#' simulation outlined in Rizopoulos (2011). For a subject last observed at time \eqn{t}, the 
#' probability that they survive until future time \eqn{u} is
#' 
#' \deqn{
#'   Pr(T_i \ge u | T \ge t; \boldsymbol{Y}_i, \boldsymbol{b}_i; \boldsymbol{\Omega}) \approx
#'   \frac{S(u|\hat{\boldsymbol{b}}_i; \boldsymbol{\Omega})}
#'   {S(t|\hat{\boldsymbol{b}}_i; \boldsymbol{\Omega})}
#' }
#' where \eqn{T_i} is the true failure time for subject \eqn{i}, \eqn{\boldsymbol{Y}_i} their
#' longitudinal measurements up to time \eqn{t}, and \eqn{S()} the survival function.
#' 
#' \eqn{\boldsymbol{\Omega}} is drawn from the multivariate normal distribution with mean
#' \eqn{\hat{\boldsymbol{\Omega}}} and its variance taken from a fitted \code{joint} object.
#' \eqn{\hat{\boldsymbol{b}}} is drawn from the \eqn{t} distribution by means of a
#' Metropolis-Hastings algorithm with \code{nsim} iterations.
#' 
#' @references 
#' 
#' Bernhardt PW, Zhang D and Wang HJ. A fast EM Algorithm for Fitting Joint Models of a Binary 
#' Response to Multiple Longitudinal Covariates Subject to Detection Limits. 
#' \emph{Computational Statistics and Data Analysis} 2015; \strong{85}; 37--53
#' 
#' Rizopoulos D. Dynamic predictions and prospective accuracy in joint models
#' for longitudinal and time-to-event data. \emph{Biometrics} 2011;
#' \strong{67}: 819–829.
#' 
#' @seealso \code{\link{ROC}} and \code{\link{plot.dynPred}}.
#' @importFrom stats median quantile
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @export
#'
#' @examples
#' \donttest{
#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#' # Focus in on id 81, who fails at around 7 years of follow-up. \code{dynPred} allows us to
#' # infer how the model believes their survival probability would've progressed (ignoring the
#' # true outcome at start time).
#' # Univariate -----------------------------------------------------------
#' long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' fit <- joint(long.formulas, surv.formula, PBC, family)
#' preds <- dynPred(PBC, id = 81, fit = fit, u = NULL, nsim = 200,
#'                  scale = 2)
#' preds
#' plot(preds)
#' # Bivariate ------------------------------------------------------------
#' # Does introduction of albumin affect conditional survival probability?
#' long.formulas <- list(
#'   serBilir ~ drug * time + I(time^2) + (1 + time + I(time^2)|id),
#'   albumin ~ drug * time + (1 + time|id)
#' )
#' fit <- joint(long.formulas, surv.formula, data = PBC, family = list("gaussian", "gaussian"))
#' bi.preds <- dynPred(PBC, id = 81, fit = fit, u = NULL, nsim = 200, 
#'                     scale = fit$coeffs$D/sqrt(fit$ModelInfo$n))
#' bi.preds
#' plot(bi.preds) # Appears to level-off dramatically; perhaps indicative of this id's albumin
#'                # levels, or acceleration in serBilir trajectory around 8.5 years.
#' }
dynPred <- function(data, id, fit, u = NULL, nsim = 200, progress = TRUE,
                    scale = NULL, df = NULL){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  
  if("factor"%in%class(data$id)) data$id <- as.numeric(as.character(data$id))
  # Check survival times
  ft <- fit$hazard[,1];tmax <- max(ft); K <- length(fit$ModelInfo$family)
  if(!is.null(u) & any(u > tmax)) stop(sprintf("Can't extrapolate beyond last failure time %.4f.\n", tmax))
  
  # Subset the required subject
  newdata <- data[data$id == id, ] # subset required subject
  
  # If u isn't supplied then arbitrarily find probability of surviving all following failure times.
  if(is.null(u)){
    last.long.time <- max(newdata$time)
    u <- c(last.long.time, ft[ft > last.long.time])
  }
  
  # Get indices for \b and \beta
  b.inds <- fit$ModelInfo$inds$C$b
  beta.inds <- fit$ModelInfo$inds$C$beta
  
  # Obtain 'denominator' dataset based on first value of vector u
  newdata2 <- newdata[newdata$time <= u[1], ]
  data.t <- prepareData(newdata2, fit, u = NULL)
  
  u <- u[-1] # Now can remove T_start

  if(nsim > 0){
    pi <- structure(matrix(NA, nrow = nsim, ncol = length(u)),
                    dimnames = list(as.character(1:nsim), paste0('u=',u)))
    MH.accept <- 0
    b.current <- shift <- data.t$bfit$par; Sigma <- solve(data.t$bfit$hessian)
    if(!is.null(scale)) Sigma <- Sigma * scale
    if(is.null(df)) df <- 4
    if(progress) pb <- utils::txtProgressBar(max = nsim, style = 3)
    for(i in 1:nsim){
      O <- Omega.draw(fit)
      b.sim <- b.mh(b.current, shift, Sigma, data.t$long, data.t$surv, O, fit, df)
      b.current <- b.sim$b.current
      MH.accept <- MH.accept + b.sim$accept
      St <- S_(data.t$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)
      for(uu in seq_along(u)){
        data.u <- prepareData(newdata, fit = fit, u = u[uu])
        pi[i, uu] <- S_(data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)/(St)
      }
      if(progress) utils::setTxtProgressBar(pb, i)
    }
    if(progress) close(pb)
    pi.df <- data.frame(
      u = u,
      mean = colMeans(pi),
      median = apply(pi, 2, median),
      lower  = apply(pi, 2, quantile, prob = 0.025),
      upper  = apply(pi, 2, quantile, prob = 0.975),
      id = id
    )
    row.names(pi.df) <- NULL
    out <- list(
      pi = pi.df,
      pi.raw = pi,
      MH.accept = MH.accept/nsim
    )
    class(out) <- 'dynPred'
    attr(out, 'type') <- 'simulated'
  }else{
    if(progress) cat("Proceeding with first-order estimate.\n")
    b <- data.t$bfit$par
    qs <- sapply(b.inds, length)
    St <- S_(data.t$surv, rep(fit$coeffs$gamma, qs), fit$coeffs$zeta, b)
    pi <- numeric(length(u))
    for(uu in seq_along(u)){
      data.u <- prepareData(newdata, fit = fit, u = u[uu])
      pi[uu] <- S_(data.u$surv, rep(fit$coeffs$gamma, qs), fit$coeffs$zeta, b)/St
    }
    pi.df <- data.frame(u = u, prob = pi, id = id)
    row.names(pi.df) <- NULL
    out <- list(pi = pi.df)
    class(out) <- 'dynPred'
    attr(out, 'type') <- 'first order'
  }
  out
} 

#' @method print dynPred
#' @keywords internal
#' @export
print.dynPred <- function(x, ...){
  if(!inherits(x, 'dynPred')) stop('x must be a "dynPred" object.')
  if(attr(x, 'type') == 'simulated'){
    cat("\nEvent-free survival probabilities based on MC simulation scheme:\n")
    print(round(x$pi[,-6], 3))
    cat(sprintf("\nMetropolis-Hastings acceptance rate %.2f%%.\n", 100 * x$MH.accept))
  }else{
    cat("\n Event-free survival probabilities based on first-order estimates:\n")
    print(round(x$pi, 3))
  }
  invisible(x)
}

#' Plot conditional survival probabilities.
#' 
#' @description Produces a simple plot of the probability the subject survives given failure
#' times along with the 95\% confidence interval if probabilities are generated from a simulation
#' scheme, and simply the first-order estimate otherwise.
#' 
#' @param x an object with class \code{dynPred}.
#' @param what should the \code{"median"} (the default) or \code{"mean"} probability be displayed?
#' if first-order predictions are used, this is ignored.
#' @param ... additional arguments (none used).
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' 
#' @importFrom graphics plot axis lines
#' @method plot dynPred
#' @seealso \code{\link{dynPred}}
#' @keywords internal
#' @export
plot.dynPred <- function(x, what = c("median", "mean"), ...){
  if(!inherits(x, 'dynPred')) stop('x must be a "dynPred" object.')
  pi <- x$pi
  id <- pi$id[1]
  if(attr(x, 'type') == 'simulated') {
    what <- match.arg(what)
    y <- pi[,what]
    sim <- TRUE
  }else{
    y <- pi[,'prob']
    sim <- FALSE
  }
  
  plot(pi$u, y, main = paste0("Event-free probability for subject ", id),
       type = 's', lwd = 1.2, col = 'black',
       xlab = expression(u), ylab = expression(Pr*"("*T[i]>=u*"|"*Y[i]*","~b[i]*";"~Omega*")"),
       ylim = 0:1, yaxt = 'n', ...)
  axis(2, at = seq(0,1,.1), labels = as.character(seq(0,1,.1)))
  if(sim) lines(pi$u, pi$lower, lwd = 1.2, lty = 3, type = 's')
  if(sim) lines(pi$u, pi$upper, lwd = 1.2, lty = 3, type = 's')
}