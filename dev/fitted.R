#' Obtain joint model fitted values
#' 
#' @description returns the fitted values from a \code{joint} object. Note that the 
#' \strong{linear predictor} for each \eqn{k=1,\dots,K} response is returned.
#'
#' @param object a joint model fit by the \code{\link{joint}} function. 
#' @param data the \emph{original} data set (i.e. that used in the \code{joint} call).
#' @param as should the fitted values be returned as a \code{"matrix"} (the default) or as a 
#' \code{"list"}? Note that \code{as="matrix"} only works for balanced responses.
#' @param ... Additional arguments (none used).
#'
#' @return A matrix (or list) with a column (or list entry) for each of the fitted linear
#' predictors.
#' @method fitted joint
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @seealso \code{\link{residuals.joint}} 
#' 
#' @export
#' @examples 
#' \donttest{
#' # Bivariate fit on PBC data -----------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin', 'platelets))
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
#' fitted(fit)
#' }
fitted.joint <- function(object, data, as = "matrix", ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(missing(data)) stop("Please provide original data in 'data' argument.")
  as <- match.arg(as, c('matrix', 'list'))
  
  
  # Re-make formulae
  M <- object$ModelInfo
  # Check and stop if responses are unbalanced and matrix is requested.
  if(as == 'matrix' & (length(M$nobs) > 1 & length(unique(M$nobs)) > 1))
    stop("Unbalanced responses, please return as = 'list' instead.")
  fs <- M$long.formulas; K <- length(fs)
  fs <- lapply(fs, parseFormula)
  
  # Data matrices
  dmats <- createDataMatrices(data, fs)
  
  # Ranefs and beta estimates
  .b <- ranef(object)
  b <- lapply(1:M$n, function(i) .b[i,,drop=F])
  b.inds <- M$inds$b
  beta <- object$coeffs$beta
  beta.inds <- M$inds$beta
  
  # Fitted value (of __linear predictor__)
  fits <- mapply(function(X, Z, b){
    lapply(1:K, function(k){
      X[[k]] %*% beta[beta.inds[[k]]] + Z[[k]] %*% b[b.inds[[k]]]
    })
  }, X = dmats$X, Z = dmats$Z, b = b, SIMPLIFY = F)
  
  # Get into one "long" column for each and return...
  out <- setNames(lapply(1:K, function(k){
    do.call(c, lapply(fits, el, k))
  }), gsub('\\s+.*$', '', M$ResponseInfo))
  if(as == 'matrix') out <- do.call(cbind, out)
  class(out) <- 'fitted.joint'
  out
}

#' Obtain joint model residuals
#' 
#' @description returns the residuals values from a \code{joint} object. The residuals
#' are of the same type as the corresponding response.
#' 
#' @param object a joint model fit by \code{\link{joint}} function.
#' @param data the \emph{original} data set (i.e. that used in the \code{joint} call).
#' @param ... Additional arguments (none used).
#' 
#' @method residuals joint
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @export
#' 
#' @seealso \code{\link{fitted.joint}}
#' 
#' @returns a named list of length \eqn{K}, with residuals produced by the joint model
#' for each of the \eqn{k=1,\dots,K} responses.
#' 
#' @examples 
#' \donttest{
#' # Trivariate fit on PBC data -----------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin', 'ascites', 'platelets'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify trivariate fit
#' long.formulas <- list(
#'   albumin ~ time*drug + (1 + time|id),
#'   platelets ~ time * drug + (1 + time|id),
#'   ascites ~ time * drug + (1|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, 
#' family = list('gaussian', 'poisson', 'binomial'))
#' residuals(fit)
#' }
residuals.joint <- function(object, data, ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(missing(data)) stop("Please provide original data in 'data' argument.")
  
  fits <- fitted(object, data, 'list') # Get fitted values
  
  M <- object$ModelInfo
  resps <- unname(sapply(M$ResponseInfo, function(x) gsub('\\s+.*$', '', x)))
  Ys <- data[,match(resps, names(data))]
  fams <- unlist(M$family); K <- length(fams)
  
  out <- setNames(lapply(1:K, function(k){
    f <- fams[k]
    switch(f,
           gaussian = Ys[,resps[k]] - fits[[resps[k]]],
           poisson = Ys[,resps[k]] - exp(fits[[resps[k]]]),
           genpois = Ys[,resps[k]] - exp(fits[[resps[k]]]),
           binomial = Ys[,resps[k]] - plogis(fits[[resps[k]]]),
           Gamma = Ys[,resps[k]] - exp(fits[[resps[k]]]))
  }), resps)
  class(out) <- "residuals.joint"
  out
}

#' Plot joint model residuals 
#' 
#' @description Produces a simple (panelled) plot of the residuals obtained by a
#' joint model (obtained by \code{\link{joint}}) of the same type as the original response.
#'
#' @param x an object with class \code{residuals.joint}.
#' @param ... additional arguments (none used).
#'
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @method plot residuals.joint
#' @importFrom graphics plot par
#' @seealso \code{\link{residuals.joint}} 
#' @export
plot.residuals.joint <- function(x, ...){
  if(!inherits(x, 'residuals.joint')) stop("Only usable with objects of class 'residuals.joint'.")
  .par <- par(no.readonly = T) # store old par
  
  K <- length(x); resps <- names(x)
  # Work out plot dimensions, set maximum of three plots per row.
  ncol <- 3
  nrow <- max(cumsum(1:K%%3 == 1))
  par(mfrow = c(nrow, ncol),
      mai = c(0.25, 0.75, 0.2, 0.25))
  
  # Plot the residuals from residuals.joint object (x).
  for(k in 1:K){
    plot(x[[k]], xaxt = 'n', main = resps[k],
         pch = 20, ylab = 'Residuals', xlab = '')
    abline(h = 0, lty = 5, col = 'red')
  }
  
  on.exit(par(.par))
}
