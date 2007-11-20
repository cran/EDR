\name{edr.R}
\alias{edr.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Eigenvectors of the effective dimension
  reduction (EDR) space }
\description{
  Computes the eigenvectors of the effective dimension
  reduction (EDR) space obtained by function \code{edr}.
}
\usage{
edr.R(B, m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{B}{ Either an object of class \code{edr} created by \code{edr} or the list component
  \code{bhat} of such an object. }
  \item{m}{ Dimension of the effective dimension reduction (EDR) space. \code{m=1} 
  corresponds to single index models, \code{m>1} specifies a multiindex model.  }
}
\value{Matrix of dimension \code{c(m,d)} containing the \code{m} eigenvectors as rows.}
}
\references{ \item{}{J. Polzehl, S. Sperlich, V. Spokoiny (2006). \emph{Estimating Generalized Principle Components} Manuscript in preparation.}}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{edr}}, \code{\link{edr.R}}, \code{\link{print.edr}}, \code{\link{summary.edr}} }
\examples{
require(EDR)
\dontrun{demo(edr_ex1)}
\dontrun{demo(edr_ex2)}
}
\keyword{ regression }% __ONLY ONE__ keyword per line
