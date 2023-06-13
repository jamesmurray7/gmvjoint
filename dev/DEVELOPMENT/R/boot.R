#' Bootstrapping a \code{joint} object
#' 
#' @description Use an existing model fit by \code{joint} along with the data object originally
#' used and obtain a mean estimate, standard errors and 95\% confidence interval using the
#' bootstrap. The original data is resampled by subject, not by observations.
#' 
#' @param fit a joint model fit by the \code{\link{joint}} function.
#' @param data the original data used to fit the above joint model.
#' @param boot.size integer, specifies the number of subjects to resample in the bootstrapping
#' approach. The default value is \code{boot.size = NULL} which defaults to the number of unique
#' subjects in the \code{joint} object.
#' @param nboot integer, specifies the number of bootstrap samples, default value is 
#' \code{nboot = 100L}.
#' @param replace logical, should sampling be done with replacement? Defaults to 
#' \code{replace = TRUE}.
#' @param progress logical, should a text progress bar showing overall progress be shown
#' and updated after each successful bootstrapped model fit? Defaults to \code{progress=TRUE}.
#' @param control a list of control arguments, with same possible arguments as shown in 
#' \code{\link{joint}}. Note that the items \code{return.dmats} \code{post.process} and
#' \code{verbose} are all set to \code{FALSE} in \code{bootjoint} in order to reduce memory
#' overheads and computation time, so these do not need to be re-specified. Instead, the user
#' could speed up this computationally intensive algorithm by changing convergence criteria
#' items e.g. \code{conv}, \code{tol.rel}, \code{tol.abs}, \code{tol.thr} in order to speed-up
#' convergence of the \code{nboot} individual bootstrapped model fits.
#' 
#' @return A list of class \code{boot.joint} which contains the MLEs from supplied \code{joint}
#' object, as well as the bootstrapped summaries and some model/computation information.
#' @seealso \code{\link{joint}} \code{\link{vcov.joint}}
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#'
#' @export
#' @examples 
#' \donttest{
#' # Bivariate fit on PBC data -----------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin', 'platelets'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify bivariate fit
#' long.formulas <- list(
#'   albumin ~ time*drug + (1 + time|id),
#'   platelets ~ time * drug + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, family = list('gaussian', 'poisson'))
#' # Set 50 bootstraps, with lower absolute tolerance and convergence of 'either'.
#' BOOT <- boot.joint(fit, PBC, nboot = 50L, control = list(tol.abs = 5e-3, conv = 'either'))
#' BOOT # Print to console via S3 method
#' }
boot.joint <- function(fit, data,
                       boot.size = NULL, nboot = 100L, replace = TRUE, progress = TRUE,
                       control = list()){
  if(!inherits(fit, 'joint')) stop("Only usable with objects of class joint.")
  if(isFALSE(fit$ModelInfo$control$post.process)) stop("Please fit parent model with post.process = TRUE.\n")
  # Get control arguments and formulae used by original joint fit.
  con <- fit$ModelInfo$control
  longs <- fit$ModelInfo$long.formulas
  disps <- fit$ModelInfo$disp.formulas
  surv <- fit$ModelInfo$surv.formula
  fam <- fit$ModelInfo$family
  if(is.null(boot.size)) boot.size <- fit$ModelInfo$n
  # Override some control arguments by default...
  con$return.dmats <- FALSE    #  don't want massive data objects.
  con$post.process <- FALSE    #  don't need to calculate hessian/posterior modes.
  con$verbose <- FALSE         #  don't want to print to console.
  # Override tolerance arguments
  conname <- names(con)
  con[(conname <- names(control))] <- control
  if(any(!names(control)%in%conname))
    warning("Double check supplied control, choice of: ",
            paste(sapply(names(fit$ModelInfo$control), sQuote), collapse = ', '), '.')
  
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
      joint(long.formulas = longs, surv.formula = surv, 
            data = resamp.data, family = fam, 
            disp.formulas = disps,
            control = con),
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
    Mean = apply(boot.coeffs, 1, mean),
    SE = apply(boot.coeffs,1,sd),
    CI = apply(boot.coeffs, 1, quantile, probs = c(.025, .975)),
    boot.time = boot.time,
    convergence.rate = sum(convs)/nboot,
    nboot = nboot,
    M = fit$ModelInfo
  ), class = 'boot.joint')
  invisible(out)
}

print.boot.joint <- function(x, digits = 3, ...){
  if(!inherits(x, 'boot.joint')) stop("Only usable with class boot.joint.\n")
  .round <- function(X) round(X, digits) # Show all to <digits> only.
  responses <- x$M$Resps       # Response names only
  K <- x$M$K                   # Number of responses
  families <- unlist(x$M$family)                     # Families only
  families <- sapply(families, neat.family.name)
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
  cat("\nLongitudinal processes:\n")
  lapply(1:K, function(k){
    cat(sprintf("\n%s (%s)\nCall:\n", responses[k], families[k]))
    cat(long.formula.to.print(x$M$long.formulas[[k]], 2), '\n')
    kresp <- responses[k]
    kinds <- which(grepl(kresp, names(x$Mean), fixed = TRUE))
    beta <- cbind(MLE = x$MLE$beta[x$M$inds$R$beta[[k]]], 
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
  print(apply(surv, 2, .round))
  
  cat("\nComputational details:\n")
  cat(sprintf("Bootstrapping took %.2f seconds\n", x$boot.time))
  cat(sprintf("%d bootstrap samples were used with %.2f%% convergence rate.", x$nboot, x$convergence.rate*100))
  cat("\n")
  invisible(x)
}
