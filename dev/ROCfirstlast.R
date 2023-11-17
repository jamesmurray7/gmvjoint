# #########################################################################
# Correcting ROC (hopefully)                                             ##
# By scheme outlined in Andrinopoulou (2021) // doi: 10.1093/ije/dyab047 ##
# #########################################################################

# Roxygen'd stuff to get into .globalenv.
#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#' long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' fit <- joint(long.formulas, surv.formula, PBC, family)
#' data <- PBC
#' Tstart <- 5; delta <- 2; control <- list()
#' progress=T; boot=T

# Create bootstrapped data
resampledata <- function(data){
  uids <- unique(data$id)
  samps <- sample(x = uids, size = length(uids), replace = TRUE)
  newData <- setNames(lapply(1:length(uids), function(i){
    newData <- data[data$id == samps[i],]
    newData$InternalKey <- i
    newData$..old.id <- newData$id
    newData$id <- i
    newData
  }),paste0('original id ', samps)) # overkill but we dont look at this anyway.
  
  as.data.frame(do.call(rbind, newData), row.names = NULL)
}





#' data(PBC)
#' PBC$serBilir <- log(PBC$serBilir)
#' long.formulas <- list(serBilir ~ drug * time + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ drug
#' family <- list('gaussian')
#' fit <- joint(long.formulas, surv.formula, PBC, family)
#' data <- PBC
#' Tstart <- 5; delta <- 2; control <- list()
#' progress=T; boot=T
corrected.ROC <- function(fit, data, Tstart, delta, control = list(), progress = TRUE,
                          nboot = 100L){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  # Parse control arguments ----
  if(!is.null(control$scale)) scale <- control$scale else scale <- 2
  if(!is.null(control$df)) df <- control$df else df <- 4
  if(!is.null(control$nsim)) nsim <- control$nsim else nsim <- 0 # 25 # Set fairly low. ## removed -> FO only atm
  sim <- nsim > 0
  
  # Ensure {survtime, status} exists via model call.
  data$survtime <- data[,fit$ModelInfo$survtime]; data$status <- data[,fit$ModelInfo$status]
  # Set out new data and remove IDs where only one longitudinal measurement is available as this causes issues in calculation 
  newdata <- data[data$survtime > Tstart, ] # subjects who are ALIVE at Tstart.
  if('factor'%in%class(newdata$id)) newdata$id <- as.numeric(as.character(newdata$id)) # this confuses tapply later on
  bad.ids <- as.numeric(names(which(with(newdata, tapply(time, id, function(x) length(unique(x)))) == 1)))
  newdata <- newdata[!newdata$id%in%bad.ids, ]
  
  # Set out candidate failure times (u)
  ft <- fit$hazard[, 1]; tmax <- max(ft)
  window <- c(Tstart + 1e-6, Tstart + 1e-6 + delta)
  if(window[2] > tmax) window[2] <- tmax
  candidate.u <- c(Tstart, ft[ft > window[1] & ft <= window[2]])
  # We solely look at their tail survival time.
  candidate.u <- candidate.u[c(1, length(candidate.u))]
  
  if(progress) pb <- utils::txtProgressBar(max = nboot, style = 3)
  for(b in 1:nboot){
    # D_(b)
    data.b <- resampledata(newdata)
    # Number of individuals alive at T_{start}
    uids <- unique(data.b$id)
    n.alive <- length(uids)
    mapped.to <- data.b[!duplicated(data.b$id), c('id', '..old.id')]
    
    probs.b <- setNames(vector('list', n.alive),
                        paste0("id ", uids)) # A bit confusing -> these _map onto_ `..old.id` in data.b
    
    pi.b <- sapply(seq_along(uids), function(i){
      ds <- dynPred(newdata, mapped.to[i,2], fit, candidate.u, progress = F,
                    scale = scale, df = df, nsim = nsim)
      out <- ds$pi
      c(orig.id = out[3], new.id = i, pi = out[2])
    })
    
  }

  
 
  
  # Loop over ids and failure times
  probs <- acceptance <- setNames(vector('list', length = length(alive.ids)),
                                  paste0('id ', as.character(unique(keys))))
  if(progress) pb <- utils::txtProgressBar(max = length(alive.ids), style = 3)
  for(i in seq_along(alive.ids)){
    ds <- dynPred(newdata, alive.ids[i], fit, u = candidate.u, progress = F, 
                  scale = scale, df = df, nsim = nsim)
    probs[[i]] <- ds$pi
    if(sim) acceptance[[i]] <- ds$MH.accept
    if(progress) utils::setTxtProgressBar(pb, i)
  }
  if(progress) close(pb)
  
  # Obtaining conditional probabilities for those alive subjects at Tstart.
  infodf <- lapply(unique(keys), function(x){
    p <- as.data.frame(probs[[paste0('id ', x)]])
    if(sim) p$mean <- p$mean else p$mean <- p$prob
    p$id <- x
    p
  })
  # pi(u|t)
  infodf <- do.call(rbind, infodf)
  pi <- with(infodf, tapply(`mean`, id, tail, 1)) # Pr(T_i^* \ge u_{max}) (i.e. survive window)
  ompi <- 1 - pi
  
  # Working out whether individuals failed/censor/survived the window
  checks <- sapply(unique(newdata$keys), function(i){
    # This key's data
    a <- newdata[newdata$keys == i,][1,,drop=F]
    survtime <- a$survtime
    status <- a$status
    # Failed or were censored otherwise
    ct.in.window <- survtime >= window[1] & survtime < window[2]
    # Did they fail?
    fail.in.window <- ct.in.window & (status == 1L)
    censor.in.window <- ct.in.window & (status == 0L)
    # Return
    c(ct.in.window = ct.in.window,
      fail.in.window = fail.in.window,
      censor.in.window = censor.in.window,
      survived.window = survtime > window[2])
  })
  # survtimes <- with(newdata, tapply(survtime, keys, unique))
  # events <- survtimes >= window[1] & survtimes <= window[2]
  # event <- status & events      # Check if they FAILED in window
  event <- checks[2,] # Specifically failed in the window
  
  # Calibration metrics
  n.window.events <- sum(event)
  # Brier score https://en.wikipedia.org/wiki/Brier_score
  #           {survived} {probs survived}   ... in the window
  BS <- mean((pi-checks[4,])^2)
  # PE.raw <- mean((checks[4,]-pi)^2) # These are the same 
  # Predictive error taking into account censoring
  PE <- mean(as.numeric(!checks[1,]) * ompi^2 + as.numeric(checks[2,]) * (0-pi)^2 + 
               as.numeric(checks[3,]) * (pi * ompi^2 + ompi * (0-pi)^2))

  # Defining threshold and calculating performance metrics.
  t <- seq(0, 1, length = 101)
  simfail <- structure(outer(pi, t, '<='),
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
    nsim = nsim,
    scale = scale,
    df = df
  )
  
  a <- AUC(out)
  out <- list(
    Tstart = Tstart, delta = delta, candidate.u = candidate.u,
    window.failures = n.window.events,
    Tstart.alive = n.alive,
    metrics = out, AUC = a, BrierScore = BS, PE = PE,
    MH.acceptance = if(sim) do.call(c, acceptance) else NULL,
    MH.acceptance.bar = if(sim) mean(do.call(c, acceptance)) else NULL,
    simulation.info = simulation.info
  )
  class(out) <- 'ROC.joint'
  out
}

# Calculate AUC from ROC table using formula for area of trapezoid.
# https://stats.stackexchange.com/a/146174
#' @keywords internal
AUC <- function(x){
  sens <- x$TPR; omspec <- x$FPR;
  height <- 0.5 * (sens[-1] + sens[-length(sens)])
  width <- -diff(omspec)
  sum(height*width)
}

#' @method print ROC.joint
#' @keywords internal
#' @export
print.ROC.joint <- function(x, ...){
  if(!inherits(x, 'ROC.joint')) stop('x must be a "ROC.joint" object.')
  q.acc <- quantile(x$MH.acceptance, c(.5,.25,.75))
  cat(sprintf("Median [IQR] M-H acceptance rate: %.3f [%.3f, %.3f]\n", q.acc[1], q.acc[2], q.acc[3]))
  cat("Diagnostic table:\n")
  print(round(x$metrics, 3))
  cat(sprintf("\nArea under curve: %.2f\n", x$AUC))
  cat(sprintf("Brier Score: %.3f\n", x$BrierScore))
  cat(sprintf("Predictive Error: %.3f\n", x$PE))
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
#' @importFrom graphics plot abline legend arrows
#' @method plot ROC.joint
#' @seealso \code{\link{dynPred}} and \code{\link{ROC}}
#' @keywords internal
#' @export
plot.ROC.joint <- function(x, legend = TRUE, show.Youden = TRUE, show.F1 = FALSE, ...){
  if(!inherits(x, 'ROC.joint')) stop('x must be a "ROC.joint" object.')
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
                  'AUC: ', round(x$AUC, 3)),
           bty = 'n', cex = .75)
  }
  invisible(x)
}

