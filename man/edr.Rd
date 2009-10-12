\name{edr}
\alias{edr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of the effective dimension reduction (EDR) space:
Structure adaptive approach for dimension reduction}
\description{
This function implements the algorithms, proposed in
M. Hristache, A. Juditsky, J. Polzehl and V. Spokoiny (2001) and 
... (2006), for estimation of the effective dimension reduction (EDR) space
in multi-index regression models \deqn{y=f(x)+\varepsilon=g(B_m^T x) + \varepsilon.} 
}
\usage{
edr(x, y, m = 2, rho0 = 1, h0 = NULL, ch = exp(0.5/max(4, (dim(x)[2]))), 
crhomin = 1, cm = 4, method = "Penalized", basis = "Quadratic", cw = NULL, 
graph = FALSE, show = 1, trace = FALSE, fx = NULL, R = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{\code{x} specifies the design matrix, dimension \code{(n,d)}  }
  \item{y}{\code{y} specifies the response, length \code{n}.}
  \item{m}{Rank of matrix M in case of \code{method="Penalized"}, not used for the other methods. }
  \item{rho0}{Initial value for the regularization parameter \eqn{\rho}. }
  \item{h0}{Initial bandwidth. }
  \item{ch}{Factor for indecreasing \eqn{h} with iterations. }
  \item{crhomin}{Factor to in(de)crease the default value of rhomin. This is just added
     to explore properties of the algorithms. Defaults to 1.  }
  \item{cm}{Factor in the definition of \eqn{\Pi_k=C_m*\rho_k^2 I_L + \hat{M}_{k-1}}. Only used if \code{method="Penalized"}. }
  \item{method}{Secifies the algoritm to use. The default \code{method="Penalized"} corresponds to the algoritm 
  proposed in ... (2006). \code{method="HJPS"} corresponds to the original algorithm from
    Hristache et.al. (2001) while \code{method="HJPS2"} specifies a modifification 
    (correction) of this algoritm. 
    }
  \item{basis}{ Specifies the set of basis functions. Options are \code{basis="Quadratic"} (default) and \code{basis="Linear"}.}
  \item{cw}{\code{cw}  another regularization parameter, secures identifiability of a minimum number of local gradient directions. Defaults to \code{1/d }. Has to be positive or \code{NULL}. }
  \item{graph}{If \code{graph==TRUE} intermediate results are plotted. }
  \item{show}{If \code{graph==TRUE} the parameter \code{show} determines the dimension of the EDR that is to be used when plotting
   intermediate results. If \code{trace=TRUE} and \code{!is.null(R)} it determines the dimension of the EDR when computing
   the risk values.}
  \item{trace}{ \code{trace=TRUE}  additional diagnostics are provided for each iteration. This includes
  current, at iteration \eqn{k}, values of the regularization parameter \eqn{\rho_k} and bandwidth \eqn{h_k},
   normalized cimmulative sums of eigenvalues of \eqn{\hat{B}} and if \code{!is.null(R)} two distances between the
   true, specified in \eqn{R} and estimated EDR.  }
  \item{fx}{ True values of \eqn{f(x)}. This is just added
     to explore properties of the algorithms and not used in the algorithms. }
  \item{R}{True matrix R. This is just added
     to explore properties of the algorithms and not used in the algorithms. }
}
\details{See reference for details.
}
\value{Object of class \code{"edr"} with components.
\item{x}{The design matrix.} 
\item{y}{The values of the response.} 
\item{bhat}{Matrix \eqn{\hat{B}} characterizing the effective dimension space. For a specified dimension \code{m} 
  \eqn{\hat{B}_m = \hat{B} O_m}, with \eqn{\hat{B}^T \hat{B}= O \Lambda O^T} being the eigenvalue decomposition of \eqn{\hat{B}^T \hat{B}}, 
  specifies the projection to the \code{m}-dimensional subspace that provides the best approximation. }
  \item{fhat}{an highly oversmoothed estimate of the values of the regression function at the design points. This is provided 
  as a backup only for the case that package \code{sm} is not installed.} 
  \item{cumlam}{Cummulative amount of information explained by the first components of \eqn{\hat{B}}.}
  \item{nmean}{Mean numbers of observations used in each iteration.}
  \item{h}{Final bandwidth}
  \item{rho}{Final value of \eqn{\rho}}
  \item{h0}{Initial bandwidth}
  \item{rho0}{Initial value of \eqn{\rho}}
  \item{cm}{The factor \code{cm}}
  \item{call}{Arguments of the call to edr}
}
\references{ M. Hristache, A. Juditsky, J. Polzehl and V. Spokoiny (2001). \emph{Structure adaptive approach for dimension reduction}, The Annals of Statistics. Vol.29, pp. 1537-1566. \\
J. Polzehl, S. Sperlich (2009). \emph{A note on structural adaptive dimension reduction},  J. Stat. Comput. Simul.. Vol. 79 (6),  pp. 805--818. 
}
\author{Joerg Polzehl, \email{polzehl@wias-berlin.de}  }
\seealso{ \code{\link{edrcv}},\code{\link{plot.edr}},   \code{\link{summary.edr}},  \code{\link{print.edr}}, \code{\link{edr.R}}}
\examples{
require(EDR)
demo(edr_ex1)
demo(edr_ex2)
}
\keyword{smooth }% at least one, from doc/KEYWORDS
\keyword{regression }% __ONLY ONE__ keyword per line
