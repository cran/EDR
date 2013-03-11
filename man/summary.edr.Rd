\name{summary.edr}
\alias{summary.edr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Summary information for objects produced by function edr. }
\description{
  The function provides information on the estimated effective dimension reduction (EDR) space.
  }
\usage{
\method{summary}{edr}(object, m = 1, R = NULL, ...)
}
\arguments{
  \item{object}{ Object of class \code{"edr"}. }
  \item{m}{ Dimension of the effective dimension reduction (EDR) space. \code{m=1} 
  corresponds to single index models, \code{m>1} specifies a multiindex model. 
  Determines the number of eigenvectors and cumulative eigenvalues to show.}
  \item{R}{ If code \code{R} specifies a matrix (dimension \code{c(k,d)}, \eqn{k>=m},
  \code{d=dim(object$x)[2]}, this matrix is interpreted as spanning the true
  EDR space. Two distances between the estimated EDR space and the space spanned 
  \code{R[1:m,]} are computed.}
  \item{...}{ Additional parameters will be ignored }
}
\details{ Provides information on the estimated effective dimension reduction (EDR) space.
The first \code{m} basis vectors and the cummulative sum of normalized eigenvalues of matrix
\code{object$bhat} are given. If \code{R} is specified the distance
\deqn{||R (I- \hat{P}_m)||/||R|| \mbox{with} \hat{P}_m = U_m^T U_m , \hat{R}_m=U_m \Lambda V^T }
and the distance specified by Li (1992) are computed.
}
\value{Returns a list with components
\item{Rhat}{(First) m eigenvectors of the estimated EDR space. }
\item{cumlam}{Cummulative sum of first m eigenvalues of \code{object$bhat}, standardized by the sum of all eigenvalues.}
\item{loss1}{If \code{R} was specified the distance \eqn{||R (I- \hat{P}_m)||/||R||} between the true and estimated
m-dimensional EDR space.}
\item{loss2}{The distance specified by Li (1992).}
}
\references{ M. Hristache, A. Juditsky, J. Polzehl and V. Spokoiny (2001). \emph{Structure adaptive approach for dimension reduction}, The Annals of Statistics. Vol.29, pp. 1537-1566. \\
J. Polzehl, S. Sperlich (2009). \emph{A note on structural adaptive dimension reduction},  J. Stat. Comput. Simul.. Vol. 79 (6),  pp. 805--818.   \\
K.-C. Li  (1992). \emph{On principal Hessian directions for data visualization and dimension reduction:
another application of Stein's lemma}, JASA, Vol. 87, pp. 1025-1039. 
}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\seealso{ \code{\link{edr}}, \code{\link{edr.R}}, \code{\link{print.edr}}, \code{\link{plot.edr}} }
\examples{
require(EDR)
\dontrun{demo(edr_ex1)}
\dontrun{demo(edr_ex2)}
}
\keyword{ smooth }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
