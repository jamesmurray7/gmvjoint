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
#' @param data the data to which the original `joint` model was fit.
#' @param id subject identifier, i.e. for which subject is the conditional survival probabilities 
#' desired?
#' @param fit a joint model fit by the  \code{joint} function. 
#' @param u a numeric \code{vector} of candidate follow-up times for which a survival probability 
#' should be calculated. Note that the first item \code{u[1]} denotes the start of the "window"
#' and is dropped from calculations. If \code{u=NULL} (the default), then the probability of 
#' surviving until the next observed failure time is returned for each longitudinal time recorded
#' for the subject \code{id}.
#' @param nsim how many Monte Carlo simulations should be carried out? Defaults to 
#' \code{nsim=200}.
#' @param progress a logical, if \code{progress=TRUE} (the default) then a progress bar displaying
#' the current percentage of simulations have been completed.
#' @param b.density character string imposing the density for the current and proposal for each
#' draw of the random effects \eqn{\boldsymbol{b}^{(l)}, l=1,\dots,\code{nsim}}. Options are 
#' \code{b.density="normal"} i.e. using the assumption that the conditional density of the random
#' effects is 
#' 
#' \deqn{\boldsymbol{b}_i^{(l)}|\boldsymbol{Y}_i(u);\boldsymbol{\Omega}^{(l)}\sim 
#'  N(\hat{\boldsymbol{b}}_i^{(l)},\hat{\boldsymbol{\Sigma}}_i^{(l)})}
#'
#' (i.e. in line with proposal by Bernhardt \emph{et al.} (2015)), or in keeping with other 
#' literature surrounding dynamic predictions (e.g. Rizopoulos (2011)) in imposing the \code{
#' b.density="t"} distribution.
#' @param scale if \code{b.density = "t"} then this numeric scales the variance-covariance 
#' parameter in the proposal distribution for the Metropolis-Hastings algorithm. Defaults to 
#' \code{scale = NULL} which doesn't scale the variance term at all.
#' @param df if \code{b.density = "t"} then this numeric denotes the degrees of freedom of the 
#' proposed \eqn{t} distribution on the random effects. \code{df=4} is suggested.
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
#' \eqn{\hat{\boldsymbol{b}}} is drawn from either the (multivariate) Normal, or \eqn{t} 
#' distribution by means of a Metropolis-Hastings algorithm with \code{nsim} iterations.
#' 
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
#' long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' fit <- joint(long.formulas, surv.formula, PBC, family, control = list(verbose=T))
#' preds <- dynPred(PBC, id = 81, fit = my.fit, u = u, nsim = 200, b.density = 'normal')
#' preds
#' plot(preds)
#' }
dynPred <- function(data, id, fit, u = NULL, nsim = 200, progress = TRUE,
                     b.density = c('normal', 't'), scale = NULL, df = NULL){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  b.density <- match.arg(b.density)
  
  # Check survival times
  ft <- fit$hazard[,1];tmax <- max(ft); K <- length(fit$ModelInfo$family)
  if(!is.null(u) & any(u > tmax)) stop("Can't extrapolate beyond last failure time.")
  
  # Subset the required subject
  newdata <- data[data$id == id, ] # subset required subject
  
  # If u isn't supplied then arbitrarily find probability of surviving until next failure time.
  #              (for each longitudinal time recorded)
  if(is.null(u)){
    u <- sapply(newdata$time, function(x) ft[which(ft > x)][1])
    if(any(u > tmax)) u <- u[!which(u > tmax)] # and ensure those after T_{max} aren't included
  }
  
  # Get indices for \b and \beta
  responsenames <- lapply(strsplit(fit$ModelInfo$ResponseInfo, '\\s\\('), el , 1)
  b.inds <- lapply(fit$ModelInfo$inds$b, function(x) x - 1)
  beta.inds <- lapply(fit$ModelInfo$inds$beta, function(x) x - 1)
  
  # Obtain 'denominator' dataset
  newdata2 <- newdata[newdata$time <= u[1], ]
  data.t <- prepareData(newdata2, id = id, fit = fit, u = NULL)
  
  u <- u[-1] # Don't want to find preds for first time T_{start}...
  
  pi <- structure(matrix(NA, nr = nsim, nc = length(u)),
                  dimnames = list(as.character(1:nsim), paste0('u=',round(u,3))))
  MH.accept <- 0
  b.current <- shift <- data.t$b$par; Sigma <- solve(data.t$b$hessian)
  if(!is.null(scale)) Sigma <- Sigma * scale
  if(b.density == 't' & is.null(df)) df <- 4
  if(progress) pb <- utils::txtProgressBar(max = nsim, style = 3)
  for(i in 1:nsim){
    O <- Omega.draw(fit)
    b.sim <- b.mh(b.current, shift, Sigma, data.t$long, data.t$surv, O, beta.inds, b.inds, fit, b.density, df)
    b.current <- b.sim$b.current
    MH.accept <- MH.accept + b.sim$accept
    St <- S_(data.t$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)
    for(uu in seq_along(u)){
      # cat('uu:', uu, '; u[uu]:', u[uu],  # uncomment for loop debugging
      #     '\nb:', b.current,'.\n')
      data.u <- prepareData(newdata, id = id, fit = fit,  u = u[uu])
      pi[i, uu] <- S_(data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)/(St)# + 1e-6)
      # cat('pi(uu):', 
      #     Surv_(data.u$surv, rep(O$gamma, sapply(b.inds, length)), O$zeta, b.current)/(St),
      #     '\n.pi[i, uu]:', pi[i,uu], '.\n')
    }
    if(progress) utils::setTxtProgressBar(pb, i)
  }
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
  out
} 

#' @method print dynPred
#' @keywords internal
#' @export
print.dynPred <- function(x, ...){
  if(!inherits(x, 'dynPred')) stop('x must be a "dynPred" object.')
  cat("\n")
  print(round(x$pi, 3))
  cat("\n")
  invisible(x)
}

#' Plot conditional survival probabilities.
#' 
#' @description Produces a simple plot of the probability the subject survives given failure
#' times along with the 95\% confidence interval.
#' 
#' @param x an object with class \code{dynPred}.
#' @param what should the \code{"median"} (the default) or \code{"mean"} probability be displayed?
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
  what <- match.arg(what)
  id <- pi$id[1]
  y <- pi[,what]
  plot(pi$u, y, main = paste0("Event-free probability for subject ", id),
       type = 's', lwd = 1.2, col = 'black',
       xlab = expression(u), ylab = expression(Pr*"("*T[i]>=u*"|"*Y[i]*","~b[i]*";"~Omega*")"),
       ylim = 0:1, yaxt = 'n')
  axis(2, at = seq(0,1,.1), labels = as.character(seq(0,1,.1)))
  lines(pi$u, pi$lower, lwd = 1.2, lty = 3, type = 's')
  lines(pi$u, pi$upper, lwd = 1.2, lty = 3, type = 's')
}