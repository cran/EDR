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
  corresponds to single index models, \code{m>1} specifies a multiindex model.}
}
\value{Matrix of dimension \code{c(m,d)} containing the \code{m} eigenvectors as rows.}
\references{ M. Hristache, A. Juditsky, J. Polzehl and V. Spokoiny (2001). \emph{Structure adaptive approach for dimension reduction}, The Annals of Statistics. Vol.29, pp. 1537-1566. \\
J. Polzehl, S. Sperlich (2009). \emph{A note on structural adaptive dimension reduction},  J. Stat. Comput. Simul.. Vol. 79 (6),  pp. 805--818. 
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{edr}}, \code{\link{edr.R}}, \code{\link{print.edr}}, \code{\link{summary.edr}} }
\examples{
require(EDR)
\dontrun{demo(edr_ex1)}
\dontrun{demo(edr_ex2)}
}
\keyword{ regression }% __ONLY ONE__ keyword per line
