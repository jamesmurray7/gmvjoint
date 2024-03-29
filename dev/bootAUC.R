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
#' @param control list of control arguments to pass to \code{\link{dynPred}} via 
#' \code{\link{ROC}}, these largely control the MC simulation scheme. Note the user
#' shouldn't specify the number of  simulations here as well, as this is covered 
#' by the argument \code{nsim}.
#'
#' @return a list of class \code{bootAUC.joint} containing relevant summary information,
#' items of interest from each call to \code{\link{ROC}} and bootstrapping settings.
#' @export
#' @seealso \code{\link{bootAUCdiff}}, \code{\link{dynPred}} and \code{\link{ROC}}.
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#'
#' @examples
#' \donttest{
#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#' long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' fit <- joint(long.formulas, surv.formula, PBC, family)
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
  cat(sprintf("Best performing probabilistic threshold by average Youden index: %.2f\n", mean(info$thresh)))
  
  invisible(x)
}

#' Compute difference in AUC between two joint models by bootstrapping.
#'
#' @description Compute summary for difference in area under the curve of two
#' competing joint model fits, given a set of data and a time-window of interest by bootstrapping.
#' 
#' @param fits a list of length two, each element containing a joint model 
#' fit by the \code{joint} function. The second element in \code{fits} should be the more complex.
#' Both elements of \code{fits} should be fit to the \strong{same} \code{data}.
#' @param data the data to which both of the \code{joint} models was fit.
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
#' @param control list of control arguments to pass to \code{\link{dynPred}} via 
#' \code{\link{ROC}}, these largely control the MC simulation scheme. Note the user
#' shouldn't specify the number of  simulations here as well, as this is covered 
#' by the argument \code{nsim}.
#'
#' @return a list of class \code{bootAUC.diff.joint} containing relevant summary information,
#' items of interest from each call to \code{\link{ROC}} and bootstrapping settings.
#' @export
#' @seealso \code{\link{bootAUC}}, \code{\link{dynPred}} and \code{\link{ROC}}.
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#'
#' @examples
#' \donttest{
#' # Compare linear vs quadratic fit to log(serBilir) in PBC data ---------
#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#' long.formulas1 <- list(serBilir ~ drug * time + (1 + time|id))
#' long.formulas2 <- list(serBilir ~ drug * (time + I(time^2)) + (1 + time + I(time^2)|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' # Fit two separate models ----------------------------------------------
#' fit1 <- joint(long.formulas1, surv.formula, PBC, family)
#' fit2 <- joint(long.formulas2, surv.formula, PBC, family)
#' # Cast to list, where more saturated/complex is second entry.
#' fits <- list(fit1, fit2)
#' # 50 (x2) first order estimates in window (8, 9] ---------------------------
#' AUC.diffs <- bootAUCdiff(fits, PBC, Tstart = 8, delta = 1, nboot = 50)
#' # print
#' AUC.diffs
#' }
bootAUCdiff <- function(fits, data, Tstart, delta, 
                        boot.size = NULL, nboot = 100, replace = TRUE, ci = 0.95,
                        nsim = 0, progress = TRUE, control = list()){
  if(!is.null(control$nsim)) 
    stop("Control number of simulations using `nsim` argument, not `control`.")
  if(!'list'%in%class(fits)) stop("'fits' must be a list of two competing fitted joint models")
  if(length(fits)!=2) stop("'fits' must contain exactly two competing models.")
  control$nsim <- nsim
  
  # If don't supply boot.size, resample the same number of ids
  # Models should be fit to SAME data.
  if(is.null(boot.size)) boot.size <- fits[[1]]$ModelInfo$n
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
  
  # Unpack list of fitted joint objects
  checks <- unlist(lapply(fits, inherits, 'joint'))
  if(sum(checks) != 2) stop("list of fits must contain 'joint' objects (failed for element) ", which(!checks))
  
  # Create empty data store lists
  info <- vector('list', nboot)
  times <- AUCs <- replicate(2, numeric(nboot), simplify = F)
  if(progress) pb <- utils::txtProgressBar(max = nboot * 2, style = 3)
  
  # Begin bootstrapping
  for(b in 1:nboot){
    bdata <- resampledata(data, boot.size, replace)
    thresh <- numeric(2)
    for(j in 1:2){
      start <- proc.time()[3]
      bROCj <- ROC(fits[[j]], bdata, Tstart, delta, control = control,
                   progress = F, boot = TRUE)
      end <- proc.time()[3]
      AUCs[[j]][b] <- bROCj$AUC
      times[[j]][b] <- end - start
      thresh[j] <- bROCj$metrics[which.max(bROCj$metrics$J)[1], "threshold"]
      if(progress) utils::setTxtProgressBar(pb, pb$getVal() + 1)
    }
    info[[b]] <- data.frame(Tstart.alive = bROCj$Tstart.alive, # This will be the same in both
                            window.failures = bROCj$window.failures,
                            thresh1 = thresh[1], thresh2 = thresh[2])
  }
  if(progress) close(pb)
  
  # Create output list
  lb <- (1 - ci)/2; ub <- 1 - lb
  AUCsumm <- lapply(AUCs, quantile, probs = c(lb, .5, ub), na.rm = TRUE)
  AUCdiff <- apply(do.call(cbind, AUCs), 1, function(x) x[2] - x[1]) # 2nd fit - 1st fit
  AUCdiffsumm <- quantile(AUCdiff, probs = c(lb, .5, ub), na.rm = TRUE)
  out <- list(
    AUCs = AUCs,
    info = info,
    ci = ci,
    AUCsumm = AUCsumm,
    AUCdiff = AUCdiff,
    AUCdiffsumm = AUCdiffsumm,
    Tstart = Tstart, delta = delta,
    boot.size, nboot = nboot, nsim = nsim
  )
  class(out) <- 'bootAUC.diff.joint'
  out
}


#' @method print bootAUC.diff.joint
#' @keywords internal
#' @export
print.bootAUC.diff.joint <- function(x, ...){
  if(!inherits(x, 'bootAUC.diff.joint')) stop('x must be a "bootAUC.diff.joint" object.')
  ci <- 100 * x$ci
  sim <- x$nsim > 0
  info <- do.call(rbind, x$info)
  cat(sprintf("Time window [%d, %d).\n", x$Tstart, x$Tstart + x$delta))
  cat(sprintf("Based on %d bootstrapped samples of the data\n\n", x$nboot))
  if(sim)
    cat("%d MC simulations used to generate survival-free probabilities.\n\n", x$nsim)
  else
    cat("First-order estimate for survival-free probabilities was used.\n\n")
  
  cat(sprintf("Model 1 median AUC [%d%% CI]:\n  %.3f [%.3f, %.3f].\n", 
              ci, x$AUCsumm[[1]][2], x$AUCsumm[[1]][1], x$AUCsumm[[1]][3]))
  
  cat(sprintf("Model 2 median AUC [%d%% CI]:\n  %.3f [%.3f, %.3f]\n\n", 
              ci, x$AUCsumm[[2]][2], x$AUCsumm[[2]][1], x$AUCsumm[[2]][3]))
  
  cat(sprintf("Median [%d%% CI] for difference in AUCs:\n  %.3f [%.3f, %.3f]\n\n",
              ci, x$AUCdiffsum[2], x$AUCdiffsum[1], x$AUCdiffsum[3]))
  
  q1 <- floor(quantile(info$Tstart.alive, probs = c(.25, .5, .75)))
  q2 <- round(quantile(info$window.failures, probs = c(.25, .5, .75)))
  cat(sprintf("Median [IQR] number of resampled subjects alive at time %d: %d [%d, %d]\n", 
              x$Tstart, q1[2], q1[1], q1[3]))
  cat(sprintf("Median [IQR] number of failures in window: %d [%d, %d].\n", q2[2], q2[1], q2[3]))
  cat(sprintf("Best-performing "))
  cat(sprintf("Best performing probabilistic threshold by average Youden index for Model 1/2: %.2f/%.2f\n", 
              mean(info$thresh1), mean(info$thresh2)))
  
  invisible(x)
}
