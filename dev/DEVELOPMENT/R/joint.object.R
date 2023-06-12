#' Fitted \code{joint} object
#'
#' @description An object returned by the \code{joint} function, with class \code{joint}
#'   a fitted joint model. Objects of this class currently have methods for: \code{logLik},
#'   \code{print}, \code{ranef}, \code{fixef}, \code{summary}, \code{AIC}, and \code{vcov}.
#'
#' @author James Murray (\email{j.murray7@@ncl.ac.uk}).
#' @seealso \code{\link{joint}}.
#' @return A list with the following components. \describe{
#'
#'  \item{\code{coeffs}}{A list containing parameter estimates: \describe{
#'  \item{\code{D}}{The variance-covariance matrix of the random effects.}
#'  \item{\code{beta}}{Vector of fixed effects for longitudinal processes.}
#'  \item{\code{sigma}}{List of dispersion parameters, families with no dispersion parameter
#'  are returned as an unnamed zero value.}
#'  \item{\code{gamma}}{Vector of association parameters.}
#'  \item{\code{zeta}}{Vector of time-invariant survival coefficients.}
#'  }}
#'  \item{\code{hazard}}{A matrix containing unique failure times \code{ft}, their hazard 
#'  contribution \code{haz} and the number of events at that failure time \code{nev}.}
#'  \item{\code{ModelInfo}}{A list containing information on the model fit: \describe{
#'  \item{\code{ResponseInfo}}{A vector containing response names with (family) reported.}
#'  \item{\code{Resps}}{A vector containing response names only.}
#'  \item{\code{family}}{A list of families fit.}
#'  \item{\code{K}}{An integer specifying the number of longitudinal sub-models.}
#'  \item{\code{Pcounts}}{A list containing informations about the number of parameters/random
#'  effects: \describe{
#'  \item{\code{P}}{A vector of length K containing the number of fixed effects for each response
#'  (in order).}
#'  \item{\code{Pd}}{A vector of length K containing the number of dispersion parameters for each
#'  response (in order) 0 denotes no parameter for that response.}
#'  \item{\code{q}}{An integer denoting the number of random effects.}
#'  \item{\code{vD}}{An integer denoting the number of unique variance-covariance parameters
#'  estimated.}
#'  }}
#'  \item{\code{long.formulas}}{A list of \code{long.formulas} (i.e. from \code{joint} call).}
#'  \item{\code{disp.formulas}}{A list of \code{disp.formulas} (i.e. from \code{joint} call).
#'  If no \code{disp.formulas} are supplied to \code{joint}, then this is populated by a list of
#'  \eqn{K} "\code{~1}". The environment is set to \code{parent.frame} in this case to avoid
#'  memory overheads in returned objects.}
#'  \item{\code{surv.formula}}{Formula object from \code{joint} call.}
#'  \item{\code{survtime}}{The name of the event time used in \code{surv.formula}.}
#'  \item{\code{status}}{The name of the event indicator used in \code{surv.formula}.}
#'  \item{\code{control}}{List of control parameters used, see \code{\link{joint}}.}
#'  \item{\code{convergence.criteria}}{List of parameters relating to the stopping rule.}
#'  \item{\code{inds}}{A list of length two, named \code{R} and \code{Cpp}, each of which contains
#'  the indices for fixed effects \eqn{\beta} for each response, or the random effects \eqn{b}
#'  for the named platform.} 
#'  \item{\code{n}}{Number of subjects.}
#'  \item{\code{nobs}}{A vector containing total number of observations for each response.}
#'  \item{\code{mi}}{A \eqn{K} x \eqn{n} matrix containing the number of observations for 
#'  subject \eqn{i} for the \eqn{k}th response.}
#'  \item{\code{nev}}{Number of events.}
#'  \item{\code{id.assign}}{A list containing the original ids of subjects in the \code{data} 
#'  supplied to \code{joint}, and the id assigned to them for use in subsequent functions.}
#'  }}
#'  \item{\code{Hessian}}{The (approximated) Hessian found at MLEs. Only returned if 
#'  control argument \code{post.process=TRUE}}.
#'  \item{\code{vcov}}{The full variance-covariance matrix between parameters. Only returned if 
#'  control argument \code{post.process=TRUE}.}
#'  \item{\code{SE}}{A named vector of approximated standard error for each estimated parameter.
#'  Only returned if control argument \code{post.process=TRUE}.}
#'  \item{\code{logLik}}{log-likelihood evaluated at parameter estimates. Only returned if control
#'  argument \code{post.process=TRUE}.}
#'  \item{\code{REs}}{The random effects, with subject-specific variance matrices attributed. If
#'  control argumnet \code{post.process=TRUE} then these are found at MLEs (i.e. are posterior
#'  estimates), otherwise they are taken from the final EM iteration.}
#'  \item{\code{elapsed.time}}{Named numeric containing breakdown of elapsed time for \code{joint}
#'  fit.}
#'  \item{\code{dmats}}{A list of data matrices on each of the longitudinal 
#'  and survival processes for each subject.}
#' }
"joint.object" <- NULL
