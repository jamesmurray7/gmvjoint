% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitted.R
\name{fitted.joint}
\alias{fitted.joint}
\title{Obtain joint model fitted values}
\usage{
\method{fitted}{joint}(object, as = "matrix", ...)
}
\arguments{
\item{object}{a joint model fit by the \code{\link{joint}} function.}

\item{as}{should the fitted values be returned as a \code{"matrix"} (the default) or as a 
\code{"list"}? Note that \code{as="matrix"} only works for balanced responses.}

\item{...}{Additional arguments (none used).}
}
\value{
A matrix (or list) with a column (or list entry) for each of the fitted linear
predictors with class \code{fitted.joint}.
}
\description{
returns the fitted values from a \code{joint} object. Note that the 
\strong{linear predictor} for each \eqn{k=1,\dots,K} response is returned.
}
\examples{
\donttest{
# Bivariate fit on PBC data -----------------------------------------
data(PBC)

# Subset data and remove NAs
PBC <- subset(PBC, select = c('id', 'survtime', 'status', 'drug', 'time',
                              'albumin', 'platelets'))
PBC <- na.omit(PBC) 

# Specify bivariate fit
long.formulas <- list(
  albumin ~ time*drug + (1 + time|id),
  platelets ~ time * drug + (1 + time|id)
)
surv.formula <- Surv(survtime, status) ~ drug

fit <- joint(long.formulas, surv.formula, PBC, family = list('gaussian', 'poisson'))
fitted(fit)
}
}
\seealso{
\code{\link{residuals.joint}}
}
\author{
James Murray (\email{j.murray7@ncl.ac.uk}).
}
