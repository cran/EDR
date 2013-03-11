\name{predict.edr}
\alias{predict.edr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Predict values of the link function for a fitted edr-object }
\description{
  The function allows to predict values of the link function in a multi-index
  model with estimated effective dimension reduction space.
}
\usage{
\method{predict}{edr}(object, xest, m = 1, h=NULL, method = "sm",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ object of class \code{edr}}
  \item{xest}{ matrix of design points where values of the link function are to be predicted }
  \item{m}{ specified dimension of the dimension reduction space. }
  \item{h}{ bandwidth, if \code{h=NULL} a bandwidth is guessed from the estimation}
  \item{method}{ \code{method="sm"} requires and uses package sm while \code{method="direct"}
                 specifies a local linear smoother. }
  \item{...}{ Additional parameters that are currently not evaluated.}
}
\details{ This function provides predictions of response values based on a multiindex model
   analyzed using function \code{edr}. It requires specification of the dimension of the
   dimension reduction space in argument \code{m}. Argument \code{method} allows to select between
   two different local smoothers for estimation of values of the link function. The use of \code{method = "sm"}
   requires package \code{sm} and is restricted to \code{m<=2}, see documentation of package sm. 
   \code{method = "direct"} uses an implementation of a local linear smoother. If parameters of the local
   linear smoother are not identified for a given bandwidth \code{h} a kernel smoother or if this is
   not identified a 1-nearest neighbor estimate is used.
The bandwidth is, in case of \code{h=NULL}, is
                 determined by generalized cross-validation for sm and as \code{x$h*x$rho}
                 in case of the local linear smoother.
}
\value{ List with components.
  \item{x }{The content of argument \code{xest}}
  \item{fhat }{Predicted values of the link function}
}
\references{M. Hristache, A. Juditsky, J. Polzehl and V. Spokoiny (2001). Structure adaptive approach for dimension reduction, \emph{The Annals of Statistics}. Vol.29, pp. 1537-1566. 

J. Polzehl, S. Sperlich (2008). A Note on Stuctural Adaptive Dimension Reduction, \emph{Journal of Statistical Computation and Simulation}, DOI: 10.1080/00949650801959699}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\note{ This function requires package \code{sm} if \code{method="sm"}.
}
\seealso{ \code{\link{edr}},\code{\link{plot.edr}},   \code{\link{summary.edr}},  \code{\link{print.edr}}, \code{\link{edr.R}}, \code{\link{edrcv}}}
\examples{
require(EDR)
\dontrun{demo(edr_ex4)}
}
\keyword{ smooth  }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
