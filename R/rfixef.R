#' Extract fixed effects from a \code{joint} object.
#' 
#' @param object a joint model fit by the \code{joint} function. 
#' @param what character string. Should the \code{"long"}itudinal process(es) be extracted,
#' or the \code{"surv"}ival ones?
#' @param ... additional arguments (none used).
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @seealso \code{\link{ranef.joint}}
#'
#' @method fixef joint
#' @return A vector containing requested fixed effects.
#' @export
#' 
#' @examples
#' # Univariate fit on PBC data -------------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify simple univariate fit
#' long.formulas <- list(
#'   albumin ~ time + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, family = list('gaussian'))
#' 
#' fixef(fit, 'long')
#' fixef(fit, 'surv')
fixef.joint <- function(object, what = c("long", 'surv'), ...){
  x <- object
  if(!inherits(x, 'joint')) stop("Only usable with objects of class 'joint'.")
  co <- x$coeffs
  what <- match.arg(what)
  if(what == 'long')
    return(c(co$beta))
  else
    return(c(co$zeta, co$gamma))
}


#' Extract random effects from a \code{joint} object.
#' 
#' @description Return the random effects \eqn{\hat{\boldsymbol{b}}} which maximises the complete
#' data log-likelihood at the MLEs \eqn{\hat{\Omega}}.
#'
#' @param object a joint model fit by the \code{joint} function. 
#' @param Var logical, should the estimated variance of the random effects at \eqn{\hat{\Omega}}
#' be returned? Defaults to \code{Var=FALSE}.
#' @param ... additional arguments (none used).
#' 
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @seealso \code{\link{fixef.joint}} \code{\link{cond.ranefs}}
#' @method ranef joint
#' @return A \code{matrix} containing required random effects effects. If \code{Var=TRUE},
#' instead a list is returned with first element the \code{matrix} of random effects and second a 
#' \code{matrix} of the variances \eqn{\hat{\Sigma}}. Note that these are \emph{posterior modes}
#' of the random effects. Conditional distribution can be found by \code{\link{cond.ranefs}}.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Univariate fit on PBC data -----------------------------------------
#' data(PBC)
#'
#' # Subset data and remove NAs
#' PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
#'                               'albumin'))
#' PBC <- na.omit(PBC) 
#' 
#' # Specify univariate fit
#' long.formulas <- list(
#'   albumin ~ time*drug + (1 + time|id)
#' )
#' surv.formula <- Surv(survtime, status) ~ drug
#' 
#' fit <- joint(long.formulas, surv.formula, PBC, family = list('gaussian'))
#' b <- ranef(fit, FALSE)
#' }
ranef.joint <- function(object, Var = FALSE, ...){
  x <- object
  if(!inherits(x, 'joint')) stop("Only usable with objects of class 'joint'.")
  if(is.null(x$REs)) stop("Rerun with post.process = TRUE.")
  RE <- x$REs
  if(Var){
    Sigma <- attr(RE, 'Var')
    attr(RE, 'Var') <- NULL
    attr(RE, 'vcov') <- NULL
    out <- list(
      b = RE,
      Sigma = Sigma
    )
    return(out)
  }else{
    attr(RE, 'Var') <- NULL
    attr(RE, 'vcov') <- NULL
    return(RE)
  }
}