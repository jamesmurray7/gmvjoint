me{gmvjoint-package}
\alias{gmvjoint-package}
\alias{gmvjoint}
\docType{package}
\title{
	\packageTitle{gmvjoint}
}
\description{
	gmvjoint allows the user to fit joint models of survival and multivariate longitudinal data. The 
	longitudinal data is specified by generalised linear mixed models (GLMMs). The joint models 
	are fit via maximum likelihood using an approximate EM algorithm first proposed by Bernhardt et
	al. (2015). The GLMMs are specified using the same syntax as for package \code{glmmTMB} Brooks et
	al. (2017). The joint models themselves are then the  flexible extensions to those in e.g.
	Wulfsohn and Tsiatis (1997). The user is able to simulate data under many different response
	types.
}

\author{ 
	  \packageAuthor{gmvjoint}
}

\references{
	Bernhardt PW, Zhang D and Wang HJ. A fast EM Algorithm for Fitting Joint Models of a Binary 
	Response to Multiple Longitudinal Covariates Subject to Detection Limits. 
	\emph{Computational Statistics and Data Analysis} 2015; \strong{85}; 37--53

	Mollie E. Brooks, Kasper Kristensen, Koen J. van Benthem, Arni Magnusson, Casper W. Berg, Anders
	Nielsen, Hans J. Skaug, Martin Maechler and Benjamin M. Bolker (2017). glmmTMB Balances Speed and
	Flexibility Among Packages for Zero-inflated Generalized Linear Mixed Modeling. 
	\emph{The R Journal}, \strong{9(2)}, 378-400.

	Murray, J and Philipson P. A fast approximate EM algorithm for joint models of survival and
	multivariate longitudinal data.\emph{Computational Statistics and Data Analysis} 2022

	Wulfsohn MS, Tsiatis AA. A joint model for survival and longitudinal data
	measured with error. \emph{Biometrics.} 1997; \strong{53(1)}, 330-339.
}
\keyword{ package }
