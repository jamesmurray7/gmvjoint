#' Obtain joint model fitted values
#' 
#' @description returns the fitted values from a \code{joint} object. Note that the 
#' \strong{linear predictor} for each \eqn{k=1,\dots,K} response is returned.
#'
#' @param object a joint model fit by the \code{\link{joint}} function. 
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
fitted.joint <- function(object, as = "matrix", ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(object$dmats)) stop("Need dmats, rerun with appropriate control arguments.")
  as <- match.arg(as, c('matrix', 'list'))
  
  M <- object$ModelInfo
  K <- M$K
  # Check and stop if responses are unbalanced and matrix is requested.
  if(as == 'matrix' & (length(M$nobs) > 1 & length(unique(M$nobs)) > 1))
    stop("Unbalanced responses, please return as = 'list' instead.")
  
  # Ranefs and beta estimates
  .b <- ranef(object)
  b <- lapply(1:M$n, function(i) .b[i,,drop=F])
  b.inds <- M$inds$R$b
  beta <- object$coeffs$beta
  beta.inds <- M$inds$R$beta
  
  dmats <- object$dmats$long

  # Fitted value (of __linear predictor__)
  fits <- mapply(function(X, Z, b){
    lapply(1:K, function(k){
      X[[k]] %*% beta[beta.inds[[k]]] + Z[[k]] %*% b[b.inds[[k]]]
    })
  }, X = dmats$X, Z = dmats$Z, b = b, SIMPLIFY = F)
  
  # Get into one "long" column for each and return...
  out <- setNames(lapply(1:K, function(k){
    do.call(c, lapply(fits, el, k))
  }), gsub('\\s+.*$', '', M$Resps))
  if(as == 'matrix') out <- do.call(cbind, out)
  class(out) <- 'fitted.joint'
  out
}

# Calculate Cox-Snell residuals.
#' @keywords internal
CoxSnellResids <- function(object){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(object$dmats)) stop("Cox Snell residuals only available with dmats in joint object")
  sv <- object$dmats$surv
  l0u <- sv$l0u; SS <- sv$SS; Fu <- sv$Fu
  b <- lapply(1:object$ModelInfo$n, function(i) object$REs[i, , drop = T]); 
  gamma <- rep(object$coeffs$gamma, sapply(object$ModelInfo$inds$R$b, length))
  zeta <- object$coeffs$zeta
  cumRisks <- mapply(function(l0u, SS, Fu, b){
    crossprod(l0u, exp(SS %*% zeta + Fu %*% (gamma * b)))
  }, l0u = l0u, SS = SS, Fu = Fu, b = b)
  cumRisks
}

#' Obtain joint model residuals
#' 
#' @description returns the Pearson residuals values from a \code{joint} object.
#' 
#' @param object a joint model fit by \code{\link{joint}} function.
#' @param what character string. Should the \code{"long"}itudinal process(es) be extracted,
#' or the \code{"surv"}ival ones?
#' @param type character. The residual type for \code{what = "long"} residuals only. Choices are
#' on the \code{"response"} scale or \code{"pearson"} residuals. Cox-Snell residuals are 
#' returend if \code{what = "surv"}.
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
#'                               'albumin', 'hepatomegaly', 'platelets'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify trivariate fit
#' long.formulas <- list(
#'   albumin ~ time*drug + (1 + time|id),
#'   platelets ~ time * drug + (1 + time|id),
#'   hepatomegaly ~ time * drug + (1|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, 
#'              family = list('gaussian', 'poisson', 'binomial'))
#' R <- residuals(fit, type = 'pearson')
#' plot(R)
#' plot(residuals(fit, what = "surv"))
#' }
residuals.joint <- function(object, what = c("longit", "surv"),
                            type = c('response', 'pearson'), ...){
  if(!inherits(object, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(object$dmats)) stop("Need dmats, rerun with appropriate control arguments.")
  type <- match.arg(type)
  what <- match.arg(what)
  
  if(what == "longit"){
    fits <- fitted(object, as = 'list') # Get fitted values
    M <- object$ModelInfo; K <- M$K
    resps <- M$Resps
    Ys <- do.call(cbind, lapply(1:K, function(k) do.call(c, lapply(object$dmats$long$Y, el, k))))
    colnames(Ys) <- resps
    Ys <- as.data.frame(Ys)
    
    fams <- unlist(M$family);
    S <- object$coeffs$sigma
    
    # Fitted values back-transformed to response scale.
    fitsk <- lapply(1:K, function(k){
      f <- fams[k]
      fitsk <-  switch(f,
                       gaussian = fits[[resps[k]]],
                       poisson = exp(fits[[resps[k]]]),
                       negbin = exp(fits[[resps[k]]]),
                       genpois = exp(fits[[resps[k]]]), 
                       binomial = plogis(fits[[resps[k]]]),
                       Gamma = exp(fits[[resps[k]]])
      )
      return(fitsk)
    })
    
    # RESPONSE residuals
    resids <- lapply(1:K, function(k){
      Ys[,k] - fitsk[[k]]
    })
    
    if(type == "pearson"){
      out <- lapply(1:K, function(k){
        f <- fams[k]; fitk <- fitsk[[k]]; res <- resids[[k]]
        r <- switch(f,
                    gaussian = res/sqrt(S[[k]]),
                    poisson = res/sqrt(fitk),
                    negbin = {
                      WW <- do.call(rbind, lapply(object$dmats$long$W, el, k))
                      phi <- exp(WW %*% S[[k]])
                      res/sqrt(fitk * (1+fitk/phi))
                    },
                    genpois = {
                      WW <- do.call(rbind, lapply(object$dmats$long$W, el, k))
                      phi <- WW %*% S[[k]]
                      res/sqrt(fitk*(1+phi)^2)
                    },
                    binomial = res/sqrt(fitk * (1 - fitk)),
                    Gamma = res/sqrt(fitk^2))
        r
      })
    }else{
      out <- resids
    }
    
    r <- setNames(lapply(1:K, function(k){
      this <- c(out[[k]]); fitk <- fitsk[[k]]
      attr(this, 'fitted') <- c(fitk)
      this
    }), resps)
    
    class(r) <- "residuals.joint"
    attr(r, 'type') <- type
    attr(r, 'what') <- what
  }else{
    r <- CoxSnellResids(object)
    class(r) <- "residuals.joint"
    attr(r, 'what') <- what
    attr(r, 'type') <- "Cox-Snell"
    attr(r, 'object') <- object
  }
  r
}

#' @method print residuals.joint
#' @keywords internal
#' @export
print.residuals.joint <- function(x, ...){
  what <- attr(x, 'what')
  type <- attr(x, 'type') 
  if(what == "longit"){
    attr(x, 'type') <- NULL; attr(x, 'class') <- NULL; attr(x, "what") <- NULL
    x <- lapply(x, function(X){
      attr(X, 'fitted') <- NULL
      c(X)
    })
    
    print(x)
    if(type == 'pearson') cat("Pearson residuals ")
    if(type == 'response') cat("Residuals ")
    cat("summary:\n")
    print(round(sapply(x, summary), 4))
    invisible(x)
  }else{
    object <- attr(x, 'object')
    attr(x, 'type') <- NULL; attr(x, 'what') <- NULL; attr(x, "object") <- NULL
    # Tis <- object$dmats$surv$Tis
    Dis <- unlist(object$dmats$ph$Delta)
    # Difference in residuals and expected
    cat("Cox-Snell residuals:\n")
    print(round(summary(x), 3))
    cat("\nMartingale residuals: \n")
    print(round(summary(Dis-x), 3))
    cat("\n")
  }
  invisible(x)
}

#' Plot joint model residuals 
#' 
#' @description Plot residuals obtained by a joint model (obtained by \code{\link{joint}}). 
#' If the \code{residuals.joint} object represents the longitudinal process, a simple (paneled)
#' plot is produced (one for each response). If the residual object contains the Cox-Snell 
#' residuals then several plots are produced (interactively): The KM estimate of survival 
#' function of said residuals and then repeated for each survival covariate in the model call
#' to \code{joint} (if requested).
#'
#' @param x an object with class \code{residuals.joint}.
#' @param strata logical, should strata (for the survival sub-model only). Defaults to 
#' \code{strata = FALSE} which produces only one plot of Cox-Snell residuals.
#' @param ... additional arguments (none used).
#'
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @method plot residuals.joint
#' @importFrom graphics plot par curve
#' @seealso \code{\link{residuals.joint}} 
#' @export
plot.residuals.joint <- function(x, strata = FALSE, ...){
  if(!inherits(x, 'residuals.joint')) stop("Only usable with objects of class 'residuals.joint'.")
  .par <- par(no.readonly = T) # store old par
  
  K <- length(x); resps <- names(x)
  what <- attr(x, 'what')
  type <- attr(x, 'type')
  if(what == 'longit'){
    # Work out plot dimensions, set maximum of three plots per row.
    if(K > 2){
      ncol <- 3
      nrow <- max(cumsum(1:K%%3 == 1))
      par(mfrow = c(nrow, ncol),
          mai = c(0.6, 0.75, 0.2, 0.25))
    }else if(K == 2){
      par(mfrow = c(1, 2),
          mai = c(0.6, 0.9, 0.2, 0.25))
    }else{
      par(mfrow = c(1,1))
    }
    # Plot the residuals from residuals.joint object (x).
    ylab <- ifelse(type == 'pearson', 'Pearson residuals', 'Residuals')
    for(k in 1:K){
      fitk <- attr(x[[k]], 'fitted')
      plot(x[[k]]~fitk, main = resps[k], cex = .75,
           pch = 20, ylab = ylab, xlab = 'Fitted')
      abline(h = 0, lty = 5, col = 'red')
    }
  }else{
    object <- attr(x, 'object')
    attr(x, 'type') <- NULL; attr(x, 'what') <- NULL; attr(x, "object") <- NULL
    sort.Ti <- sort(object$dmats$surv$Tis)
    # KM estimate
    sf.null <- survfit(Surv(x, unlist(object$dmats$ph$Delta)) ~ 1)
    # Plot Survival function of residuals, and then for each strata level
    plot(sf.null, mark.time = F, 
         xlab = "Cox-Snell residuals", ylab = "Survival probability",
         main = "Survival function of Cox-Snell residuals")
    curve(exp(-x), from = 0, to = max(sort.Ti), add = T, lwd = 2, col = 'steelblue')
    
    # Stratified by call to `joint`.
    if(strata){
      xx <- readline("Press Enter for residual plots by strata: ")
      num.strata <- length(object$dmats$ph$invar.surv.names)
      for(s in 1:num.strata){
        this.strata <- object$dmats$ph$invar.surv.names[s]
        ff <- as.formula(paste0("Surv(x, status) ~ ", this.strata))
        sf.fitd <- survfit(ff, data = object$dmats$ph$survdata)
        # Work out what the data labels are
        if(object$ModelInfo$control$center.ph){
          s.strata.names <- c(names(sf.fitd$strata))
          if(length(s.strata.names) > 5) next # Don't bother (i.e. s is continuous)
          to.match <- round(as.numeric(gsub(paste0(this.strata, '='), '', s.strata.names)), 4)
          raw <- object$dmats$ph$ph$x[,this.strata]
          rs <- cbind(raw = raw, scaled = round(c(scale(raw, scale = F)), 4))
          rs <- rs[!duplicated.matrix(rs),]
          for.legend <- paste0(this.strata, " = ", rs[match(to.match, rs[, "scaled"]), "raw"])
        }else{
          s.strata.names <- c(names(sf.fitd$strata))
          if(length(s.strata.names) > 5) next # Don't bother (i.e. s is continuous)
          for.legend <- c(names(sf.fitd$strata))
        }
        # Make the plot
        plot(sf.fitd, lty = 1:2, 
             xlab = "Cox-Snell residuals", ylab = "Survival probability",
             main = "Survival function of Cox-Snell residuals by strata")
        curve(exp(-x), from = 0, to = max(sort.Ti), add = T, lwd = 2, col = 'steelblue')
        legend('topright', lty = c(1:num.strata, 1), col = c(rep('black', 2), "steelblue"),
               lwd = c(rep(1,num.strata), 2), bty = 'n',
               legend = c(for.legend, "exp(-x)"))
        # Make next plot (if possible)
        if(s < num.strata) xx <- readline("Press Enter for next strata.")
      }
    }
  }
  on.exit(par(.par))
}


