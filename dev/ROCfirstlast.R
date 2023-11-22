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

fix.hazard <- function(fit.orig, fit.boot, mapping, w, v){
  b.orig <- lapply(1:fit.orig$ModelInfo$n, function(i) fit.orig$REs[i,])
  S.orig <- lapply(1:fit.orig$ModelInfo$n, function(i) vech2mat(attr(fit.orig$REs, 'vcov')[i,],
                                                                fit.orig$ModelInfo$Pcounts$q))
  gamma.rep <- rep(fit.boot$coeffs$gamma, sapply(fit.orig$ModelInfo$inds$Cpp$b, length))
  l0.hat.boot <- c(lambda_hat(b.orig, fit.orig$dmats$surv$Fu, fit.orig$dmats$surv$SS, S.orig, gamma.rep,
                              fit.boot$coeffs$zeta, fit.orig$dmats$ph$nev, w, v))
  l0u <- lapply(fit.orig$dmats$surv$l0u, function(x) l0.hat.boot[1:length(x)])
  l0i <- lapply(1:fit.orig$ModelInfo$n, function(i){
    if(fit.orig$dmats$ph$Delta[[i]] == 1L)
      return(l0.hat.boot[match(fit.orig$dmats$surv$Tis[i], fit.orig$dmats$surv$ft)])
    else
      return(0)
  })
  
  # We need b.hat and Sigma.hat to calculate lambda based on Omega at 
  # bootstrapped parameter estimates. Chance that subject i won't have been 
  # randomly sampled in generating D_(b), in this case we need to work it out
  # in same fashion as approx EM.
  b.sig <- lapply(1:fit.orig$ModelInfo$n, function(i){
    # Check if `id` i is in the bootstrapped sample
    check <- mapping[,'..old.id'] == i
    if(any(check)){ # If they are, get b.hat and Sigma.hat for the subject from M_(b)
      # Potentially multiple assignments, so just take the first
      wh <- which(check)
      boot.id <- unname(mapping[wh[1], 'id'])
      b.hat <- fit.boot$REs[boot.id,]
      Sigma.hat <- vech2mat(attr(fit.boot$REs, 'vcov')[boot.id,], fit.boot$ModelInfo$Pcounts$q)
      return(list(b = b.hat, S = Sigma.hat))
    }else{ # Otherwise we need to calculate b.hat and Sigma.hat for them
      b.update <- optim(rep(0, fit.orig$ModelInfo$Pcounts$q), 
                        joint_density, joint_density_ddb,
                        Y = fit.orig$dmats$long$Y[[i]], X = fit.orig$dmats$long$X[[i]], 
                        Z = fit.orig$dmats$long$Z[[i]], W = fit.orig$dmats$long$W[[i]], 
                        beta = fit.boot$coeffs$beta, D = fit.boot$coeffs$D, sigma = fit.boot$coeffs$sigma, 
                        family = fit.orig$ModelInfo$family, 
                        Delta = fit.orig$dmats$ph$Delta[[i]], S = fit.orig$dmats$surv$S[[i]], 
                        Fi = fit.orig$dmats$surv$Fi[[i]], l0i = l0i[[i]], 
                        SS = fit.orig$dmats$surv$SS[[i]], Fu = fit.orig$dmats$surv$Fu[[i]], haz = l0u[[i]], 
                        gamma_rep = gamma.rep, zeta = fit.boot$coeffs$zeta,
                        beta_inds = fit.orig$ModelInfo$inds$Cpp$beta, b_inds = fit.orig$ModelInfo$inds$Cpp$b, 
                        K = fit.orig$ModelInfo$K, method = 'BFGS', hessian = TRUE)
      return(list(b = b.update$par, S = solve(b.update$hessian)))
    }
  })
  b.hats <- lapply(b.sig, '[[', 1)
  S.hats <- lapply(b.sig, '[[', 2)
  lambda.update <- c(lambda_update(b.hats, fit.orig$dmats$surv$Fu, fit.orig$dmats$surv$SS,
                                   S.hats, fit.orig$dmats$surv$surv.times, gamma.rep, fit.boot$coeffs$zeta,
                                   fit.orig$dmats$surv$nev, w, v))
  l0u <- lapply(fit.orig$dmats$surv$l0u, function(ll){
    lambda.update[1:length(ll)]
  })
  l0i <- lambda.update[match(fit.orig$dmats$surv$Ti, fit.orig$dmats$surv$ft)] 
  l0i[is.na(l0i)] <- 0
  l0i <- as.list(l0i)
  return(list(l0u = l0u, l0i = l0i, b.hats = b.hats, S.hats = S.hats))
}

# Getting the prognostic performance measures; below is for my sanity because getting confused.
# This could be coded/done a lot more efficiently: Just getting something that works for now...
#' @param fit (either) the original model, or the model fit to bootstrapped sample.
#' @param data (either) the original data, or the bootstrapped data.
get.measures <- function(fit, data, Tstart, delta, nsim){
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
  pi <- lapply(seq_along(uids), function(i){
    ds <- dynPred(newdata, uids[i], fit, candidate.u, progress = F,
                  scale = scale, df = df, nsim = nsim)
    out <- unname(ds$pi)
    data.frame(orig.id = out[3], new.id = i, pi = out[2])
  })
  pi.df <- do.call(rbind, pi)
  pi <- pi.df$pi
  ompi <- 1 - pi
    
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
  
  # w, v
  GH <- statmod::gauss.quad.prob(fit$ModelInfo$control$gh.nodes, "normal",
                                 0, fit$ModelInfo$control$gh.sigma)
  w <- GH$w; v <- GH$n
  
  # Get apparent AUC and PE on M
  M.measures <- get.measures(fit, data, Tstart, delta, nsim = nsim)
  
  if(progress) pb <- utils::txtProgressBar(max = nboot, style = 3)
  fit.copy <- fit
  opts <- matrix(NA, nrow = nboot, ncol = 2)
  for(b in 1:nboot){
    # D_(b)
    data.b <- resampledata(data)
    mapping <- simplify2array(data.b[!duplicated(data.b$id), c('id', '..old.id')])
    
    # Refit the model and get M_(b) on D_(b)
    M.b <- refit.model(fit, data.b)
    
    # Obtain bootstrapped performance measures
    M.b.measures <- get.measures(M.b, data.b, Tstart, delta, nsim)
    
    # Obtain bootstrapped performance measures on original data
    fit.copy$coeffs <- M.b$coeffs
    fit.copy$Hessian <- M.b$Hessian
    fit.copy$vcov <- M.b$vcov
    b.step <- fix.hazard(fit, M.b, mapping, w, v)
    fit.copy$dmats$surv$l0i <- b.step$l0i
    fit.copy$dmats$surv$l0u <- b.step$l0u
    b.REs <- do.call(rbind, b.step$b.hats)
    attr(b.REs, 'vcov') <- do.call(rbind, lapply(b.step$S.hats, vech))
    fit.copy$REs <- b.REs
    M.bstar.measures <- get.measures(fit.copy, data, Tstart, delta, nsim)
    
    # Calculate the optimism in AUC and PE
    opts[b,] <- unlist(M.b.measures) - unlist(M.bstar.measures)
    
    if(progress) utils::setTxtProgressBar(pb, b)
  }
  
  corrected.AUCs <- M.measures$auc - opts[,1]
  corrected.PEs <- M.measures$PE + opts[,2]
  
  out <- structure(list(AUC.c = corrected.AUCs, PE.c = corrected.PEs,
                        M.auc = M.measures$auc, M.PE = M.measures$PE,
                        nboot = nboot, Tstart = Tstart, delta = delta),
                   class = 'corrected.joint')
  out
}

print.corrected.joint <- function(x, ...){
  stopifnot(inherits(x, 'corrected.joint'))
  aa <- sprintf("Time Window w: [%.1f, %.1f],", x$Tstart, x$Tstart + x$delta)
  bb <- sprintf("Based on %d bootstrap model fits and discrimination measures.", x$nboot)
  auc <- x$AUC.c; pe <- x$PE.c
  cc <- ''
  if(any(auc>1)){
    cc <- sprintf("There was %d corrected AUC value >1, which has been removed", sum(auc > 1))
    auc <- auc[which(auc < 1)]
    pe <- pe[which(auc < 1)]
  }
  qn.auc <- quantile(auc, c(.25,.5,.75))
  qn.pe <- quantile(pe, c(.25,.5,.75))
  
  cat(aa, '\n')
  cat(bb, '\n')
  if(nchar(cc)) cat(cc, "\n")
  cat("\nAUC ----\n")
  cat(sprintf("'Point estimate': %.3f\n", x$M.auc))
  cat(sprintf("Median [IQR]: %.3f [%.3f, %.3f]\n\n", qn.auc[2], qn.auc[1], qn.auc[3]))
  cat("PE ----\n")
  cat(sprintf("'Point estimate': %.3f\n", x$M.PE))
  cat(sprintf("Median [IQR]: %.3f [%.3f, %.3f]\n\n", qn.pe[2], qn.pe[1], qn.pe[3]))
  invisible(x)
}
