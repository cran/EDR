\name{plot.edr}
\alias{plot.edr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Plot results produced by function edr. }
\description{
  Illustrate the fitted model within the estimated one or two-dimensional effective dimension
  reduction (EDR) space. 
  }
\usage{
plot.edr(x, m = 1, ylab = "Y", title = "", sm = require(sm), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Object of class \code{"edr"}. }
  \item{m}{ Dimension of the effective dimension reduction (EDR) space. \code{m=1} 
  corresponds to single index models, \code{m>1} specifies a multiindex model. 
  Values of \code{m>2} are currently not allowed.}
  \item{ylab}{ Label for the response variable. }
  \item{title}{ Optional title. }
  \item{sm}{ If \code{sm==TRUE} nonparametric regression is performed
  within the m-dimensional EDR using function \code{sm.regression} from
  package \code{sm}. If \code{sm==FALSE} or \code{require(sm)==FALSE} 
  oversmoothed fitted values calculated within function \code{edr} are used
  for visualisation. }
  \item{...}{ Additional parameters will be passed to \code{plot} in case of \code{m=1} and
  to \code{image} if \code{m=2}.}
}
\value{Returns \code{invisible{NULL}}.
}
\references{ \item{}{J. Polzehl, S. Sperlich, V. Spokoiny (2006). \emph{Estimating Generalized Principle Components} Manuscript in preparation.}}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{edr}}, \code{\link{edr.R}}, \code{\link{print.edr}}, \code{\link{summary.edr}} }
\examples{
require(EDR)
\dontrun{demo(edr_ex1)}
\dontrun{demo(edr_ex2)}
}
\keyword{ hplot }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
