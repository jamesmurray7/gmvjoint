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
#' predictors with class \code{fitted.joint}.
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
#' fitted(fit)
#' }
fitted.joint <- function(object, data = NULL, as = "matrix", ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(object$dmats) && is.null(data))
    stop("Must provide original 'data' if 'joint' object doesn't contain data matrices.")
  as <- match.arg(as, c('matrix', 'list'))
  
  M <- object$ModelInfo
  K <- length(M$family)
  # Check and stop if responses are unbalanced and matrix is requested.
  if(as == 'matrix' & (length(M$nobs) > 1 & length(unique(M$nobs)) > 1))
    stop("Unbalanced responses, please return as = 'list' instead.")
  
  # Ranefs and beta estimates
  .b <- ranef(object)
  b <- lapply(1:M$n, function(i) .b[i,,drop=F])
  b.inds <- M$inds$b
  beta <- object$coeffs$beta
  beta.inds <- M$inds$beta
  
  if(is.null(object$dmats)){ # If return.dmats is FALSE...
    # Re-make formulae
    fs <- M$long.formulas; K <- length(fs)
    fs <- lapply(fs, parseFormula)
    
    dmats <- createDataMatrices(data, fs)
  }else{
    dmats <- object$dmats$long
  }
  
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
#' @description returns the Pearson residuals values from a \code{joint} object.
#' 
#' @param object a joint model fit by \code{\link{joint}} function.
#' @param data the \emph{original} data set (i.e. that used in the \code{joint} call).
#' @param type character. The residual type
#' @param ... Additional arguments (none used).
#' 
#' @method residuals joint
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @export
#' 
#' @seealso \code{\link{fitted.joint}} \code{\link{plot.residuals.joint}}
#' 
#' @returns a named list of length \eqn{K} of class \code{residuals.joint} containing
#' residuals produced by the joint model for each of the \eqn{k=1,\dots,K} responses, 
#' along with the fitted values as an attribute.
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
#' R <- residuals(fit, type = 'pearson')
#' plot(R)
#' }
residuals.joint <- function(object, data = NULL, type = c('response', 'pearson'), ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(object$dmats) && is.null(data))
    stop("Must provide original 'data' if 'joint' object doesn't contain data matrices.")
  type <- match.arg(type)
  
  fits <- fitted(object, data = data, as = 'list') # Get fitted values
  
  M <- object$ModelInfo; K <- length(M$family)
  resps <- unname(sapply(M$ResponseInfo, function(x) gsub('\\s+.*$', '', x)))
  if(!is.null(data)){
    Ys <- data[,match(resps, names(data))]
  }else{
    # Re-construct data from dmats
    Ys <- do.call(cbind, lapply(1:K, function(k) do.call(rbind, lapply(object$dmats$long$Y, el, k))))
    colnames(Ys) <- resps
    Ys <- as.data.frame(Ys)
  }
  
  fams <- unlist(M$family);
  S <- unlist(object$coeffs$sigma)
  
  out <- setNames(lapply(1:K, function(k){
    f <- fams[k]
    fitsk <-  switch(f,
                     gaussian = fits[[resps[k]]],
                     poisson = exp(fits[[resps[k]]]),
                     genpois = exp(fits[[resps[k]]]), 
                     binomial = plogis(fits[[resps[k]]]),
                     Gamma = exp(fits[[resps[k]]])
    )
    res <- Ys[,resps[k]] - fitsk
    r <- switch(type,
                response = res,
                pearson = {
                  switch(f, 
                         gaussian = res/sqrt(S[k]),
                         poisson = res/sqrt(fitsk),
                         genpois = res/sqrt(fitsk*(1+S[k])^2),
                         binomial = res/sqrt(fitsk * (1 - fitsk)),
                         Gamma = res/sqrt(fitsk^2)
                         )
                }
    )
    list(fitted = fitsk, residuals = r)
  }), resps)
  r <- lapply(out, el, 2)
  class(r) <- "residuals.joint"
  attr(r, 'fitted') <- lapply(out, el, 1)
  attr(r, 'type') <- type
  r
}

#' @method print residuals.joint
#' @keywords internal
#' @export
print.residuals.joint <- function(x, ...){
  type <- attr(x, 'type') 
  attr(x, 'type') <- NULL; attr(x, 'fitted') <- NULL
  attr(x, 'class') <- NULL
  print(x)
  if(type == 'pearson') cat("Pearson residuals ")
  if(type == 'response') cat("Residuals ")
  cat("summary:\n")
  print(round(sapply(x, summary), 4))
  invisible(x)
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
  type <- attr(x, 'type')
  # Work out plot dimensions, set maximum of three plots per row.
  if(K > 2){
    ncol <- 3
    nrow <- max(cumsum(1:K%%3 == 1))
    par(mfrow = c(nrow, ncol),
        mai = c(0.6, 0.75, 0.2, 0.25))
  }else{
    par(mfrow = c(1, 2),
        mai = c(0.6, 0.9, 0.2, 0.25))
  }
  # Plot the residuals from residuals.joint object (x).
  ylab <- ifelse(type=='pearson', 'Pearson residuals', 'Residuals')
  for(k in 1:K){
    fitk <- attr(x, 'fitted')[[resps[k]]]
    plot(x[[k]]~fitk, main = resps[k], cex = .75,
         pch = 20, ylab = ylab, xlab = 'Fitted')
    abline(h = 0, lty = 5, col = 'red')
  }
  
  on.exit(par(.par))
}
