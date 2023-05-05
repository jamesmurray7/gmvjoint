bootjoint <- function(fit, data,
                      boot.size = NULL, nboot = 100L, replace = TRUE, progress = TRUE){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class joint.")
  # Get control arguments and formulae used by original joint fit.
  con <- fit$ModelInfo$control
  longs <- fit$ModelInfo$long.formulas
  surv <- fit$ModelInfo$surv.formulas
  fam <- fit$ModelInfo$family
  if(is.null(boot.size)) boot.size <- fit$ModelInfo$n
  
  # Function to resample data ----
  resampledata <- function(data, size, replace){
    uids <- unique(data$id)
    samps <- sample(x = uids, size = size, replace = replace)
    newData <- setNames(lapply(1:size, function(i){
      newData <- data[data$id == samps[i],]
      newData$InternalKey <- i
      newData$..old.id <- newData$id
      newData$id <- i
      newData
    }),paste0('original id ', samps)) # overkill but we dont look at this anyway.
    
    as.data.frame(do.call(rbind, newData))
  }
  
  # The bootstrapping wrapper ----
  if(progress) pb <- utils::txtProgressBar(max = nboot, style = 3)
  # Initialise stores
  boot.start <- proc.time()[3]
  boot.fits <- vector("list", nboot)
  for(b in 1:nboot){
    resamp.data <- resampledata(data, boot.size, replace)
    fit.b <- tryCatch(
      joint(longs, surv, resamp.data, fam, post.process = FALSE,
            control = c(con, return.dmats = FALSE)),
      error = function(e) NULL
    )
    if(!is.null(fit.b)) boot.fits[[b]] <- list(Omega = fit.b$coeffs, time = fit.b$elapsed.time)
    utils::setTxtProgressBar(pb, b)
  }
  boot.time <- proc.time()[3] - boot.start
  convs <- !sapply(boot.fits, is.null)
  
  if(progress) close(pb)
  
  boot.coeffs <- sapply(boot.fits, function(x){
    O <- x$Omega
    vD <- vech(O$D); beta <- O$beta; sigma <- unlist(O$sigma)[unlist(O$sigma)!=0L]
    gamma <- O$gamma; zeta <- O$zeta
    setNames(c(vD, beta, sigma, gamma, zeta), names(fit$SE))
  }) 
  boot.coeffs <- boot.coeffs[!rowSums(boot.coeffs)==0,]
  
  # Construct the output
  out <- structure(list(
    MLE = fit$coeffs,
    Mean = apply(boot.coeffs,1,mean),
    SE = apply(boot.coeffs,1,sd),
    CI = apply(boot.coeffs, 1, quantile, probs = c(.025, .975)),
    boot.time = boot.time,
    convergence.rate = sum(convs)/nboot,
    nboot = nboot,
    M = fit$ModelInfo
  ),class = 'boot.joint')
  invisible(out)
}

print.boot.joint <- function(x, digits = 3, ...){
  if(!inherits(x, 'boot.joint')) stop("Only usable with class boot.joint.\n")
  .round <- function(X) round(X, digits) # Show all to <digits> only.
  responses <- lapply(sapply(M$ResponseInfo,       # Response names only
                             strsplit, '\\s\\('), el, 1)
  K <- length(M$ResponseInfo)                      # Number of responses
  families <- unlist(M$family)                     # Families only
  cat("MLE estimates, with Bootstrap mean, SE and 95% CI estimates.\n\n")
  
  # Covariance
  cat("vech(D):\n")
  vD.inds <- which(grepl("^D\\[", names(x$SE)))
  vD.nms <- names(x$SE)[vD.inds]
  vD <- cbind(MLE = vech(x$MLE$D), Mean = x$Mean[vD.inds], SE = x$SE[vD.inds],
              `2.5%` = x$CI[1,vD.inds], `97.5%` = x$CI[2,vD.inds])
  w <- max(nchar(rownames(vD))) + 2 + max(nchar(format(.round(vD[,1]), nsmall = digits)))
  cat(rep("", w), 'Bootstrap ----\n')
  print(apply(vD,2,.round))
  
  # Longitudinal
  K <- length(x$M$family)
  cat("\nLongitudinal processes:\n")
  lapply(1:K, function(k){
    cat(sprintf("\n%s (%s)\nCall:\n", responses[[k]], families[k]))
    cat(deparse(x$M$long.formulas[[k]]), '\n')
    kresp <- responses[[k]]
    kinds <- which(grepl(kresp, names(x$Mean)))
    beta <- cbind(MLE = x$MLE$beta[x$M$inds$beta[[k]]], 
                  Mean = x$Mean[kinds], SE = x$SE[kinds],
                  `2.5%` = x$CI[1,kinds], `97.5%` = x$CI[2,kinds])
    w <- max(nchar(rownames(beta))) + 2 + max(nchar(format(.round(beta[,1]), nsmall = digits)))
    cat(rep("", w), "Boostrap----\n")
    print(apply(beta,2,.round))
  })
  
  # Survival
  cat("\nEvent-time sub-model:\n")
  cat(sprintf("Call: %s\n", deparse(x$M$surv.formula)))
  gz.MLE <- c(x$MLE$gamma, x$MLE$zeta)
  gz.nm <- names(gz.MLE)
  gz.inds <- match(gz.nm, names(x$Mean))  
  surv <- cbind(MLE = gz.MLE, 
                Mean = x$Mean[gz.inds], SE = x$SE[gz.inds],
                `2.5%` = x$CI[1,gz.inds], `97.5%` = x$CI[2,gz.inds])
  w <- max(nchar(rownames(surv))) + 2 + max(nchar(format(.round(surv[,1]), nsmall = digits)))
  cat(rep("", w), "Boostrap----\n")
  print(apply(surv,2,.round))
  
  cat("\nComputational details:\n")
  cat(sprintf("Bootstrapping took %.2f seconds\n", x$boot.time))
  cat(sprintf("%d bootstrap samples were used with %.2f convergence rate.", x$nboot, x$convergence.rate))
  cat("\n")
  invisible(x)
}