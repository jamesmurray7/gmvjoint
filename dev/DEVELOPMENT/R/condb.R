#' Obtain conditional distribution of the random effects
#'
#' @description Obtain the conditional distribution of the random effects of a \code{joint} model
#' fit. This is achieved by a Metropolis scheme. Approximate normality across random effects is
#' expected, and could be useful in diagnosing potential issues surrounding model fits.
#'
#' @param fit a joint model fit by the \code{joint} function. 
#' @param burnin Number of burn-in iterations to discard.
#' @param N Number of MC iterations to carry out post burn-in.
#' @param tune Tuning parameter, problem-specific.
#'
#' @return A list of class \code{cond.b.joint} containing: \describe{
#'  \item{walks}{A list of length \code{n} containing the history of \eqn{b_i} post burn-in.}
#'  \item{acceptance}{A numeric vector containing the acceptance rate for each sampled subject.}
#'  \item{M}{The ModelInfo list from \code{joint}. Used by S3 methods for class 
#'           \code{cond.b.joint}.}
#'  \item{bhats}{Posterior estimates at MLEs for the random effects. Same as \code{ranef(joint)}.}
#'  \item{q}{Dimension of random effects.}
#'  \item{K}{Number of responses.} 
#'  \item{qnames}{The names of the random effects as determined by call to \code{joint}.}
#'  \item{burnin}{The amount of burn-in used.}
#'  \item{N}{Number of MC iterations.}
#'  \item{tune}{tuning parameter used}
#'  \item{nobs}{The number of observations for each subject for each response.}
#'  \item{elapsed.time}{Time taken for \code{cond.ranefs} to complete.}
#' }
#' @export
#' @seealso \code{\link{ranef.joint}}
#'
#' @examples
#' \donttest{
#' dat <- simData()$data
#' long.formulas <- list(Y.1 ~ time + cont + bin + (1 + time|id), 
#'                       Y.2 ~ time + cont + bin + (1 + time|id))
#' surv.formula <- Surv(survtime, status) ~ bin
#' fit <- joint(long.formulas, surv.formula, dat, list("gaussian","gaussian"))
#' cond.b <- cond.ranefs(fit, burnin = 50L, N = 1000, tune = 2)
#' cond.b
#' plot(cond.b)
#' }
cond.ranefs <- function(fit, burnin = 500L, N = 3500L, tune = 2.){
  if(!inherits(fit, "joint")) stop("Only usable with object of class 'joint'.")
  if(is.null(fit$dmats)) stop("Need dmats, rerun with appropriate control arguments.")
  # Unpack dmats
  M <- fit$ModelInfo
  dm <- fit$dmats$long
  sv <- fit$dmats$surv
  surv <- fit$dmats$ph
  q <- sv$q; K <- M$K
  
  # Unpack parameter estimates
  Omega <- fit$coeffs
  gamma.rep <- rep(Omega$gamma, sapply(M$inds$R$b, length))
  
  # Model matrices
  l0i <- sv$l0i; l0u <- sv$l0u
  Del <- surv$Delta
  Fu <- sv$Fu; Fi <- sv$Fi
  S <- sv$S; SS <- sv$SS
  X <- dm$X; Y <- dm$Y; Z <- dm$Z
  # Other
  beta.inds <- M$inds$Cpp$beta
  b.inds <- M$inds$Cpp$b
  ff <- M$family
  nobs <- sapply(lapply(Y, el), length)
  
  # Inits + metropolis scheme inputs
  b <- lapply(1:M$n, function(i) c(fit$REs[i,,drop=F]))
  Sigmas <- lapply(1:M$n, function(i) attr(fit$REs, 'vcov')[i,])
  Sigmas <- lapply(Sigmas, vech2mat, q)
  iters <- burnin + N
  out <- vector('list', M$n); accepts <- numeric(M$n)
  pb <- utils::txtProgressBar(max = M$n, style = 3)
  start.time <- proc.time()[3]
  for(i in 1:M$n){
    this <- metropolis(b[[i]], Omega, Y[[i]], X[[i]], Z[[i]], ff, 
                       Del[[i]], S[[i]], Fi[[i]], l0i[[i]], SS[[i]], 
                       Fu[[i]], l0u[[i]], gamma.rep, beta.inds, b.inds, 
                       K, q, burnin, N, Sigmas[[i]], tune)
    utils::setTxtProgressBar(pb, i)
    out[[i]] <- t(this$walks)
    accepts[i] <- this$AcceptanceRate
  }
  end.time <- proc.time()[3]
  close(pb)
  out <- list(walks = out, acceptance = accepts, M = M, bhats = do.call(rbind, b),
              q = q, K = K, qnames = colnames(fit$coeffs$D),
              burnin = burnin, N = N, tune = tune, nobs = nobs,
              elapsed.time = end.time-start.time)
  class(out) <- 'cond.b.joint'
  out
}

#' @keywords internal
#' @method print cond.b.joint
print.cond.b.joint <- function(x, ...){
  if(!inherits(x, 'cond.b.joint')) stop("x must be a 'cond.b.joint' object.")
  M <- x$M
  # Data information
  cat("Conditional distribution f(b_i|Y_i,T_i,Delta_i;Omega)\n")
  cat("------------------------------\n")
  cat(sprintf("Number of subjects: %d\n", M$n))
  cat(sprintf("Number of events: %d (%.2f%%)\n", M$nev, 100 * M$nev/M$n))
  cat(sprintf("Number of responses: %d, dimension of random effects: %d\n", x$K, x$q))
  nobs <- x$nobs
  qnobs <- quantile(nobs)
  if(qnobs[3] == qnobs[2] && qnobs[3] == qnobs[4])
    cat(sprintf("Median profile length: %d\n", qnobs[3]))
  else
    cat(sprintf("Median [IQR] profile length: %d [%d, %d]\n", qnobs[3], qnobs[2], qnobs[4]))
  
  # MC Information
  cat("\n")
  accepts <- x$acceptance
  walks <- do.call(rbind, x$walks)
  colnames(walks) <- x$qnames
  cat("Metropolis Algorithm information:\n")
  cat(sprintf("Burn-in: %d, MC iterations post burn-in: %d\n", x$burnin, x$N))
  cat(sprintf("Median acceptance rate: %.2f%%\n", median(accepts) * 100))
  
  # Print summary of the things 
  cat("\n")
  cat("Summary: \n")
  summ <- apply(walks, 2, summary)
  summ <- rbind(summ[1:4,], Emp.SD = apply(walks, 2, sd), summ[5:6,])
  print(round(summ, 4))
  
  invisible(x)
}

#' @keywords internal
#' @method plot cond.b.joint
plot.cond.b.joint <- function(x, D = NULL, nrow = NULL, ncol = NULL, title = NULL, ...){
  if(!inherits(x, 'cond.b.joint')) stop("x must be a 'cond.b.joint' object.")
  .par <- par(no.readonly = T)
  walks <- do.call(rbind, x$walks)
  num.plots <- x$q # Determine number of plots
  if(is.null(nrow) & !is.null(ncol) | !is.null(nrow) & is.null(ncol)) 
    stop("Define either both 'nrow' and 'ncol', or leave both unspecified.")
  if(is.null(nrow) & is.null(ncol)){ # Plot as if intercept and slope by default.
    nrow <- max(cumsum(unlist(x$M$ind$b)%%2))
    ncol <- 2
  }
  par(mfrow = c(nrow, ncol))
  for(j in 1:num.plots){
    dj <- density(walks[,j])
    if(!is.null(D)){
      seqx <- seq(min(dj$x), max(dj$x), length.out = 1e3)
      dn <- dnorm(seqx, mean = 0, sd = sqrt(D[j,j]))
      ylims <- c(0, max(max(dj$y), max(dn)) + .1)
    }else{
      ylims <- c(min(dj$y), max(dj$y) + 0.1)
    }
    plot(dj, xlab = '', main = x$qnames[j],
         ylim = ylims)
    if(!is.null(D))
      lines(seqx, dn, lty = 3, col = 'red2')
  }
  if(!is.null(title))
    mtext(title, side = 3, line = -2, outer = TRUE)
  on.exit(par(.par))
}


