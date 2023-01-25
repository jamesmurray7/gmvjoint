#' Obtain bootstrapped area under the curve measures for a joint model.
#'
#' @description Compute summary for the area under the curve of a joint model fit
#' given a set of data and a time-window of interest by bootstrapping.
#' 
#' @param fit a joint model fit by the  \code{joint} function.
#' @param data the data to which the original \code{joint} model was fit.
#' @param Tstart The start of the time window of interest, \code{Tstart} denotes the time
#' point up to which longitudinal process(es) is used in calculation of survival probabilities.
#' @param delta scalar denoting the length of time interval to check for failure times.
#' @param boot.size the number of samples to take from the data. Default value is \code{NULL}
#' which uses the number of distinct \code{ids} in the \code{data}.
#' @param nboot the number of bootstrap replicates upon which the summaries are formed.
#' @param replace logical, can \code{id}s be randomly sampled more than once? Defaults to
#' \code{replace = TRUE}.
#' @param ci confidence interval to be returned on AUC estimates, defaults to 
#' \code{ci = 0.95} which produces a 95\% confidence interval.
#' @param nsim number of simulations to use in each call to \code{\link{dynPred}},
#' defaults to \code{nsim=0} which results in first-order estimates being used.
#' @param progress logical, should a progress bar showing the percentage of completed
#'  bootstrapped estimates? Defauts to \code{progress=TRUE}.
#' @param DP.control list of control arguments to pass to \code{\link{dynPred}} via 
#' \code{\link{ROC}}, these largely control the MC simulation scheme. Note the user
#' shouldn't specify the number of  simulations here as well, as this is covered 
#' by the argument \code{nsim}.
#'
#' @return a list of class \code{bootAUC.joint} containing relevant summary information,
#' items of interest from each call to \code{\link{ROC}} and bootstrapping settings.
#' @export
#' @seealso \code{\link{dynPred}} and \code{\link{ROC}}.
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#'
#' @examples
#' \donttest{
#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#' long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' fit <- joint(long.formulas, surv.formula, PBC, family, control = list(verbose=F))
#' # 100 first order estimates in window (8, 9]
#' AUCs <- bootAUC(fit, PBC, Tstart = 8, delta = 1)
#' AUCs
#' }
bootAUC <- function(fit, data, Tstart, delta, 
                    boot.size = NULL, nboot = 100, replace = TRUE, ci = 0.95,
                    nsim = 0, progress = TRUE, control = list()){
  if(!is.null(control$nsim)) 
    stop("Control number of simulations using `nsim` argument, not `control`.")
  control$nsim <- nsim
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  
  # If don't supply boot.size, resample the same number of ids
  if(is.null(boot.size)) boot.size <- fit$ModelInfo$n
  resampledata <- function(data, size, replace){
    uids <- unique(data$id)
    samps <- sample(x = uids, size = size, replace = replace)
    newData <- setNames(lapply(1:size, function(i){
      newData <- data[data$id == samps[i],]
      newData$InternalKey <- i
      newData
    }),paste0('original id ', samps))
    
    as.data.frame(do.call(rbind, newData))
  }
  
  if(progress) pb <- utils::txtProgressBar(max = nboot, style = 3)
  info <- vector('list', nboot)
  AUCs <- times <- numeric(nboot)
  for(b in 1:nboot){
    bdata <- resampledata(data, boot.size, replace)
    start <- proc.time()[3]
    bROC <- ROC(fit, bdata, Tstart, delta, control = control,
                progress = F, boot = TRUE)
    end <- proc.time()[3]
    AUCs[b] <- bROC$AUC
    info[[b]] <- data.frame(Tstart.alive = bROC$Tstart.alive,
                            window.failures = bROC$window.failures,
                            thresh = bROC$metrics[which.max(bROC$metrics$J)[1],'threshold'])
    times[b] <- end - start
    if(progress) utils::setTxtProgressBar(pb, b)
  }
  if(progress) close(pb)
  
  # Create output list
  lb <- (1 - ci)/2; ub <- 1 - lb
  qu <- quantile(AUCs, probs = c(lb, 0.5, ub)) # lower, med, upper
  
  out <- list(
    AUCs = AUCs,
    info = info,
    ci = ci,
    qu = qu,
    mean = mean(AUCs), sd = sd(AUCs),
    Tstart = Tstart, delta = delta,
    boot.size, nboot = nboot, nsim = nsim
  )
  class(out) <- 'bootAUC.joint'
  out
}


#' @method print bootAUC.joint
#' @keywords internal
#' @export
print.bootAUC.joint <- function(x, ...){
  if(!inherits(x, 'bootAUC.joint')) stop('x must be a "bootAUC.joint" object.')
  ci <- 100 * x$ci
  sim <- x$nsim > 0
  info <- do.call(rbind, x$info)
  cat(sprintf("Time window [%d, %d).\n", x$Tstart, x$Tstart + x$delta))
  cat(sprintf("Based on %d bootstrapped samples of the data\n\n", x$nboot))
  if(sim)
    cat("%d MC simulations used to generate survival-free probabilities.\n", x$nsim)
  else
    cat("First-order estimate for survival-free probabilities was used.\n")
  
  cat(sprintf("Median [%d%% confidence interval] for AUC:\n", ci))
  cat(sprintf("%.3f [%.3f, %.3f]\n\n", x$qu[2], x$qu[1], x$qu[3]))
  
  q1 <- floor(quantile(info$Tstart.alive, probs = c(.25, .5, .75)))
  q2 <- round(quantile(info$window.failures, probs = c(.25, .5, .75)))
  cat(sprintf("Median [IQR] number of resampled subjects alive at time %d: %d [%d, %d]\n", 
              x$Tstart, q1[2], q1[1], q1[3]))
              
  cat(sprintf("Median [IQR] number of failures in window: %d [%d, %d].\n", q2[2], q2[1], q2[3]))
  cat(sprintf("Best performing probabilistic threshold by average Youden index: %.2f", mean(info$thresh)))
  
  invisible(x)
}