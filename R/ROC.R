#' Receiver Operator Characteristics (ROC) for a \code{joint} model.
#' 
#' @description Using longitudinal information available up to a time, establish diagnostic
#' capabilities (ROC, AUC and Brier score) of a fitted joint model.
#'
#' @param fit a joint model fit by the  \code{joint} function.
#' @param data the data to which the original \code{joint} model was fit.
#' @param Tstart The start of the time window of interest, \code{Tstart} denotes the time
#' point up to which longitudinal process(es) is used in calculation of survival probabilities.
#' @param delta scalar denoting the length of time interval to check for failure times.
#' @param control list of control arguments to be passed to \code{\link{dynPred}}, which
#' acts as the main workhorse function for \code{ROC}. Takes default arguments of 
#' \code{\link{dynPred}} if not supplied.
#'
#' @return A list of class \code{ROC.joint} consisting of: \describe{
#'   \item{\code{Tstart}}{numeric denoting the start of the time window of interest; all dynamic
#'   predictions generated used longitudinal information up-to time \eqn{T_{\text{start}}}.}
#'   \item{\code{delta}}{scalar which denotes length of interval to check, such that the window
#'   is defined by \eqn{[T_{\text{start}}, T_{\text{start}}, + \delta]}.}
#'   \item{\code{candidate.u}}{candidate vector of failure times to calculate dynamic probability
#'    of surviving for each subject alive in \code{data} at time \eqn{T_{\text{start}}}.}
#'   \item{\code{window.failures}}{numeric denoting the number of observed failures in
#'   \eqn{[T_{\text{start}}, T_{\text{start}}, + \delta]}.}
#'   \item{\code{Tstart.alive}}{numeric denoting the risk set at \code{Tstart}.}
#'   \item{\code{metrics}}{a \code{data.frame} containing probabilistic \code{thresholds} with:
#'   \code{TP} true positives; \code{FN} false negatives; \code{FP} false positives;
#'   \code{TN} true negatives; \code{TPR} true positive rate (sensitivity); \code{FPR} false
#'   positive rate (1-specificity); \code{Acc} accuracy; \code{PPV} positive predictive value
#'   (precision); \code{NPV} negative predictive value; \code{F1s} F1 score and \code{J} Youden's
#'   J statistic.}
#'   \item{AUC}{the area under the curve.}
#'   \item{BrierScore}{calculated Brier score (for each subject) along with attributed summary.}
#'   \item{MH.acceptance.bar}{mean acceptance of M-H scheme across all subjects.}
#'   \item{simulation.info}{list containing information about call to \code{dynPred}.}
#' }
#' @export
#' @seealso \code{\link{dynPred}} and \code{\link{plot.ROC.joint}}
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
#' roc <- ROC(fit, PBC, Tstart = 8, delta = 2, control = list(nsim = 25))
#' roc
#' plot(roc)
#' }
ROC <- function(fit, data, Tstart, delta, control = list()){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  
  # Parse control arguments ----
  if(!is.null(control$b.density)) b.density <- control$b.density else b.density <- 'normal'
  b.density <- match.arg(b.density, c('normal', 't'))
  if(!is.null(control$scale)) scale <- control$scale else scale <- NULL
  if(b.density == 't' & is.null(scale)){
    message('Scale not supplied for t distribution, defaulting to * 2')
    scale <- 2
  }
  if(!is.null(control$df)) df <- control$df else df <- NULL
  if(b.density == 't' & is.null(df)){
    message('df not supplied for t distribution, defaulting to df = 4')
    df <- 4
  }
  if(!is.null(control$nsim)) nsim <- control$nsim else nsim <- 25 # Set fairly low.
  
  # Set out new data and remove IDs where only one longitudinal measurement is available as this causes issues in calculation 
  newdata <- data[data$survtime > Tstart, ] # subjects who are ALIVE at Tstart.
  if('factor'%in%class(newdata$id)) newdata$id <- as.numeric(as.character(newdata$id)) # this confuses tapply later on
  bad.ids <- as.numeric(names(which(with(newdata, tapply(time, id, function(x) length(unique(x)))) == 1)))
  newdata <- newdata[!newdata$id%in%bad.ids, ]
  alive.ids <- unique(newdata$id)
  n.alive <- length(alive.ids)
  
  # Set out candidate failure times (u)
  ft <- fit$hazard[, 1]; tmax <- max(ft)
  window <- c(Tstart + 1e-6, Tstart + 1e-6 + delta)
  if(window[2] > tmax) window[2] <- tmax
  candidate.u <- c(Tstart, ft[ft > window[1] & ft <= window[2]])
  
  # Loop over ids and failure times
  probs <- acceptance <- setNames(vector('list', length = length(alive.ids)),
                                  paste0('id ', alive.ids))
  pb <- utils::txtProgressBar(max = length(alive.ids), style = 3)
  for(i in seq_along(alive.ids)){
    ds <- dynPred(newdata, alive.ids[i], fit, u = candidate.u, progress = F, 
                   b.density = b.density, scale = scale, df = df, nsim = nsim)
    probs[[i]] <- ds$pi#do.call(rbind, ds$pi)
    acceptance[[i]] <- ds$MH.accept
    utils::setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # Obtaining conditional probabilities for those alive subjects at Tstart.
  infodf <- lapply(alive.ids, function(x){
    p <- as.data.frame(probs[[paste0('id ', x)]])
    p$mean <- p$mean   
    p$id <- x
    p
  })
  # pi(u|t)
  infodf <- do.call(rbind, infodf)
  pi <- with(infodf, tapply(`mean`, id, min))
  
  # Working out whether individuals failed in the window
  survtimes <- with(newdata, tapply(survtime, id, unique))
  events <- survtimes >= window[1] & survtimes <= window[2]
  status <- as.logical(with(newdata, tapply(status, id, unique))) # Check if they failed
  event <- status & events # Check if they failed in window.
  
  n.window.events <- sum(event)
  BS <- ((1 - pi) - sapply(event, as.numeric))^2   # Brier score
  attr(BS, 'summary') <- summary(BS)
  
  # Defining threshold and calculating performance metrics.
  t <- seq(0, 1, length = 101)
  simfail <- structure(outer(pi, t, '<'),
                       dimnames = list(names(pi) , paste0('t: ', t)))
  
  TP <- colSums(c(event) * simfail)        # True positives
  FN <- sum(event) - TP                    # False negatives
  FP <- colSums(c(!event) * simfail)       # False positives
  TN <- sum(!event) - FP                   # True negatives
  TPR <- TP/(TP + FN)                      # True positive rate (sensitivity)
  FPR <- FP/(FP + TN)                      # False positive rate (1 - specificity)
  Acc <- (TP + TN) / (TP + TN + FP + FN)   # Accuracy
  PPV <- TP/(TP + FP + 1e-6)               # Positive predictive value (precision)
  NPV <- TN/(TN + FN + 1e-6)               # Negative predictive value
  F1s <- 2*(PPV* TPR) / (PPV + TPR + 1e-6) # F1 score
  J <- TPR + TN/(TN + FP) - 1              # Youden's J statistic
  
  # Sanity checks -- mainly for internal use.
  if(!all.equal(TPR, TP/sum(event))) stop('Something wrong: TP + FN != sum(event)')
  if(!all.equal(FP / (FP + TN)  , FP/(sum(!event)))) stop('Something wrong: FP + TN != sum(!event)')
  
  # Making a nice dataframe to report
  out <- data.frame(threshold = t,
                    TP = TP, TN = TN, FP = FP, FN = FN,
                    TPR = TPR, FPR = FPR, PPV = PPV, NPV = NPV,
                    Acc = Acc, F1 = F1s, J = J)
  row.names(out) <- NULL
  
  # Flip table so if multiple thresholds have same TPR/FPR then we take the largest threshold
  out <- out[order(out$threshold, decreasing = T), ]
  # Remove duplicated TPR/FPR
  out <- out[!duplicated(out[, c('TPR', 'FPR')]),]
  
  simulation.info <- list(
    b.density = b.density,
    nsim = nsim,
    scale = scale,
    df = df
  )
  
  a <- AUC(out)
  out <- list(
    Tstart = Tstart, delta = delta, candidate.u = candidate.u,
    window.failures = n.window.events,
    Tstart.alive = n.alive,
    metrics = out, AUC = a, BrierScore = BS,
    MH.acceptance.bar = mean(do.call(c, acceptance)),
    simulation.info = simulation.info
  )
  class(out) <- 'ROC.joint'
  out
}

# Calculate AUC from ROC table.
#' @keywords internal
AUC <- function(x){
  TPR <- rev(x$TPR); FPR <- rev(x$FPR);
  sum(0.5 * diff(FPR) * (TPR[-1] + TPR[-length(TPR)]), na.rm = TRUE)
}

#' @method print ROC.joint
#' @keywords internal
#' @export
print.ROC.joint <- function(x, ...){
  if(!inherits(x, 'ROC.joint')) stop('x must be a "ROC.joint" object.')
  cat("Diagnostic table:\n")
  print(round(x$metrics, 3))
  cat(sprintf("\nArea under curve: %.2f\n", x$AUC))
  cat(sprintf("Median Brier Score: %.3f\n", unname(attr(x$BrierScore, 'summary')['Median'])))
  invisible(x)
}


#' Plot receiver operator characteristics.
#' 
#' @description Produces a simple plot showing the true positive rate (sensitivity) against
#' the false positive rate (1-specificy) for a dynamic prediction routine on a \code{joint} model
#' along a specified time interval.
#' 
#' @param x an object with class \code{ROC.joint}.
#' @param legend should a legend displaying the number in risk set; number of failures in interval;
#' area under the ROC curve and median Brier score be added to the bottom-right corner of the ROC 
#' plot? Default is \code{legend = TRUE}.
#' @param show.Youden should a line be drawn showing optimal cut-point using Youden's J statistic?
#' Defaults to \code{show.Youden = TRUE}.
#' @param show.F1 should a line be drawn showing optimal cut-point using the F-score?
#' Defaults to \code{show.F1 = FALSE}. Note that this measure comes under heavy criticism and is
#' included for completeness' sake.
#' @param ... additional arguments (none used).
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' 
#' @importFrom graphics plot abline legend arrorws
#' @method plot ROC.joint
#' @seealso \code{\link{dynPred}} and \code{\link{ROC}}
#' @keywords internal
#' @export
plot.ROC.joint <- function(x, legend = TRUE, show.Youden = TRUE, show.F1 = FALSE, ...){
  TPR <- x$metrics$TPR; FPR <- x$metrics$FPR
  plot(FPR, TPR,
       xlab = '1 - Specificity', ylab = 'Sensitivity',
       main = paste0('ROC curve for time interval (', x$Tstart, ', ', x$Tstart + x$delta, ']'),
       type = 'l')
  abline(0, 1, lty = 3)
  if(show.Youden){
    Ms <- x$metrics
    maxJ <- max(Ms$J); ind <- which.max(Ms$J)
    arrows(x0 = Ms[ind, 'FPR'], x1 = Ms[ind, 'FPR'],
           y0 = Ms[ind, 'FPR'], y1 = Ms[ind, 'TPR'],
           length = 0, lty = 5)
  }
  if(show.F1){
    Ms <- x$metrics
    maxF1 <- max(Ms$F1); ind <- which.max(Ms$F1)
    arrows(x0 = Ms[ind, 'FPR'], x1 = Ms[ind, 'FPR'],
           y0 = Ms[ind, 'FPR'], y1 = Ms[ind, 'TPR'],
           length = 0, lty = 3)
  }
  if(legend){
    legend('bottomright', 
           paste0(x$Tstart.alive, ' at risk; ', x$window.failures, ' failures in interval.\n',
                  'AUC: ', round(x$AUC, 3),
                  ', Brier score: ', round(unname(attr(x$BrierScore, 'summary')['Median']), 4), '.'),
           bty = 'n', cex = .75)
  }
  invisible(x)
}

