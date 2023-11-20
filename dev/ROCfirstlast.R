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
#' progress=T; boot=T; nboot = 10L
#' 


# Setting out some helper functions ---------------------------------------
# Creating bootstrapped data
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

# Refitting the model on the bootstrapped data
refit.model <- function(fit, D.b){ # M: The original model (for info); D.b bootstrapped data
  # Extract formulae and families
  long <- fit$ModelInfo$long.formulas
  surv <- fit$ModelInfo$surv.formula
  fam <- fit$ModelInfo$family
  disp <- fit$ModelInfo$disp.formulas
  
  # Fit with tryCatch wrapper
  tryCatch(
    joint(long.formulas = long,
          surv.formula = surv,
          data = D.b,
          family = fam, disp.formulas = disp, # For sake of comp. speed, increase tolerance and start close to MLEs
          control = list(conv = 'either', tol.abs = 5e-3, tol.rel = 1e-2,
                         inits = fit$coeffs)),
    error = function(e) NULL
  )
}

# dynPred -> additional arguments to account for M.b -> original data.
dynPred2 <- function(data, orig.id, boot.id, fit, u = NULL, nsim = 200, progress = TRUE,
                     scale = NULL, df = NULL){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  
  if("factor"%in%class(data$id)) data$id <- as.numeric(as.character(data$id))
  # Check survival times
  ft <- fit$hazard[,1];tmax <- max(ft); K <- length(fit$ModelInfo$family)
  if(!is.null(u) & any(u > tmax)) stop(sprintf("Can't extrapolate beyond last failure time %.4f.\n", tmax))
  
  # Subset the required subject
  newdata <- data[data$id == orig.id, ] # subset required subject
  newdata$id <- boot.id                 # They will appear in the `dmat`s by _this_ `id`
  
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
      id = orig.id
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
    pi.df <- data.frame(u = u, prob = pi, orig.id = orig.id, boot.id = boot.id)
    row.names(pi.df) <- NULL
    out <- list(pi = pi.df)
    class(out) <- 'dynPred'
    attr(out, 'type') <- 'first order'
  }
  out
}

# Getting the prognostic performance measures; below is for my sanity because getting confused.
# This could be coded/done a lot more efficiently: Just getting something that works for now...
#' @param fit (either) the original model, or the model fit to bootstrapped sample.
#' @param data (either) the original data, or the bootstrapped data.
#' @param mapping data.frame with the ids as they appear in the bootstrapped data vs the original.
get.measures <- function(fit, data, mapping = NULL, Tstart, delta, nsim){
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
  
  # Number of individuals alive at T_{start}
  uids <- unique(newdata$id)
  n.alive <- length(uids)
  
  # Getting pi and (1 - pi), dependent on whether mapping has occurred.
  if(is.null(mapping)){ # If no mapping has occurred (i.e. orig -> orig or boot -> boot)
    pi <- lapply(seq_along(uids), function(i){
      ds <- dynPred(newdata, uids[i], fit, candidate.u, progress = F,
                    scale = scale, df = df, nsim = nsim)
      out <- unname(ds$pi)
      data.frame(orig.id = out[3], new.id = i, pi = out[2])
    })
    pi.df <- do.call(rbind, pi)
    pi <- pi.df$pi
    ompi <- 1 - pi
  }else{ # If mapping has occurred (i.e. we're working out boot -> orig)
    # Was the original id in the original sample at T_{start}?
    mapping2 <- mapping[which(mapping[, 2]%in%uids),]
    boot.id <- mapping2[, 1]  # Their `id`s as they appear in `fit` dmats
    orig.id <- mapping2[, 2]  # Their original `id`s
    pi <- lapply(seq_along(orig.id), function(i){
      ds <- dynPred2(newdata, orig.id[i], boot.id[i], fit, candidate.u, progress = F,
                    scale = scale, df = df, nsim = nsim)
      out <- unname(ds$pi)
      data.frame(orig.id = orig.id[i], new.id = boot.id[i], pi = out[2])
    })
    pi.df <- do.call(rbind, pi)
    pi <- pi.df$pi
    ompi <- 1 - pi
  }
  # Working out whether individuals failed/censor/survived the window
  checks <- sapply(pi.df$orig.id, function(i){
    # This key's data
    a <- newdata[newdata$id == i,][1,,drop=F]
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
  
  event <- checks[2,] # Specifically failed in the window
  
  PE <- mean(as.numeric(!checks[1,]) * ompi^2 + as.numeric(checks[2,]) * (0-pi)^2 + 
               as.numeric(checks[3,]) * (pi * ompi^2 + ompi * (0-pi)^2))
  
  t <- seq(0, 1, length = 101)
  simfail <- structure(outer(pi, t, '<='),
                       dimnames = list(names(pi) , paste0('t: ', t)))
  
  TP <- colSums(c(event) * simfail)        # True positives
  FN <- sum(event) - TP                    # False negatives
  FP <- colSums(c(!event) * simfail)       # False positives
  TN <- sum(!event) - FP                   # True negatives
  TPR <- TP/(TP + FN)                      # True positive rate (sensitivity)
  FPR <- FP/(FP + TN)                      # False positive rate (1 - specificity)
  
  df <- data.frame(threshold = t, TPR = TPR, FPR = FPR)
  # Flip table so if multiple thresholds have same TPR/FPR then we take the largest threshold
  df <- df[order(df$threshold, decreasing = T), ]
  # Remove duplicated TPR/FPR
  df <- df[!duplicated(df[, c('TPR', 'FPR')]),]
  
  .auc1 <- 0.5 * (df$TPR[-1] + df$TPR[-length(df$TPR)])
  .auc2 <- -diff(df$FPR)
  # Return AUC and PE
  return(list(auc = sum(.auc1*.auc2), PE = PE))
}

corrected.ROC <- function(fit, data, Tstart, delta, control = list(), progress = TRUE,
                          nboot = 100L){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class 'joint'.")
  # Parse control arguments ----
  if(!is.null(control$scale)) scale <- control$scale else scale <- 2
  if(!is.null(control$df)) df <- control$df else df <- 4
  if(!is.null(control$nsim)) nsim <- control$nsim else nsim <- 0 # 25 # Set fairly low. ## removed -> FO only atm
  sim <- nsim > 0
  
  # Get apparent AUC and PE 
  M.measures <- get.measures(fit, data, NULL, Tstart, delta)
  
  if(progress) pb <- utils::txtProgressBar(max = nboot, style = 3)
  for(b in 1:nboot){
    # D_(b)
    data.b <- resampledata(data)
    mapping <- simplify2array(data.b[!duplicated(data.b$id), c('id', '..old.id')])
    
    # Refit the model
    M.b <- refit.model(fit, data.b)
    
    # Obtain bootstrapped performance measures
    M.b.measures <- get.measures(M.b, data.b, NULL, Tstart, delta)
    
    # Obtain bootstrapped performance measures on original data
    M.bstar.measures <- get.measures(M.b, data, mapping, Tstart, delta)
    
    # Calculate the optimism
    O.b <- 
    
  }
  
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

