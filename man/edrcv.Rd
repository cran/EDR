\name{edrcv}
\alias{edrcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Risk assessment by Cross-Validation }
\description{Tis function, additionally to estimating the
effective dimension reduction space (EDR), see also function \code{\link{edr}}, estimates the Mean Squared Error of Prediction (MSEP) and the Mean Absolute Error of Prediction (MAEP) when using the estimated EDR by Cross-Validation. Estimates of the regression function are produced using function \code{sm.regression}  from package \code{sm}.
}
\usage{
edrcv(x, y, m = 2, rho0 = 1, h0 = NULL, ch = exp(0.5/max(4, (dim(x)[2]))), crhomin = 1, 
      cm = 4, method = "Penalized",fit = "sm" , basis = "Quadratic", cw = NULL, 
      graph = FALSE, show = 1, trace = FALSE, seed = 1, cvsize = 1, m0 = min(m, 2), 
      hsm = NULL)
}
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
  \item{fit}{Specifies the method for estimating and predicting values of the link function. This can either be \code{fit="sm"}
    specifying use of the sm package or \code{fit="direct"} specifying the use of a local linear smoother. In case of \code{m0>2} 
    \code{fit="direct"} is used due to restrictions in the sm package. }
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
  \item{seed}{ Seed for generating random groups for CV }
  \item{cvsize}{ Groupsize k in leave-k-out CV  }
  \item{m0}{ Dimension of the dimension reduction space to use when fitting 
   the data. Should be either 1 or 2.  }
  \item{hsm}{ If \code{is.null(hsm)} the bandwidth used by \code{sm.regression} for smoothing within the EDR is chosen by cross-validation within \code{sm.regression} when needed. Alternatively a grid of bandwidths may be specified. In that case a bandwidth for  \code{sm.regression} is chosen from the grid that minimizes the extimated mean absolute error of prediction.  }
}
\details{
This function performs a leave-k-out cross-validation to estimate the risk
in terms of Mean Squared Error of Prediction (MSEP) and Mean Absolute Error of Prediction (MAEP) when using function \code{edr} to estimate an 
effective dimension reduction space of dimension \code{m0} and using this estimated space to predict values of the response. Smoothing within the dimension reduction space is performed using  the  function \code{sm.regression}  from package \code{sm}. The bandwidth for  \code{sm.regression} is 
chosen by Cross-Validation. 
}
\value{
Object of class \code{"edr"} with components.
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
  \item{call}{Arguments of the call to edrcv}
  \item{cvres}{Residuals from cross-validation.}
  \item{cvmseofh}{Estimates of MSEP for bandwidths \code{hsm}}
  \item{cvmaeofh}{Estimates of MAEP for bandwidths \code{hsm}}
  \item{cvmse}{Estimate of MSEP}
  \item{cvmae}{Estimate of MAEP}
  \item{hsm}{Set of bandwidths specified for use with \code{sm.regression}}
  \item{hsmopt}{Bandwidth selected for use with \code{sm.regression} if \code{hsm} was specified.}
  
}
\references{M. Hristache, A. Juditsky, J. Polzehl and V. Spokoiny (2001). Structure adaptive approach for dimension reduction, \emph{The Annals of Statistics}. Vol.29, pp. 1537-1566. 

J. Polzehl, S. Sperlich (2008). A Note on Stuctural Adaptive Dimension Reduction, \emph{Journal of Statistical Computation and Simulation}, DOI: 10.1080/00949650801959699}
\author{ Joerg Polzehl, \email{polzehl@wias-berlin.de} }
\note{ This function requires package \code{sm} if \code{fit="sm"}.
}
\seealso{ \code{\link{edr}},\code{\link{plot.edr}},   \code{\link{summary.edr}},  \code{\link{print.edr}}, \code{\link{edr.R}}, \code{\link{predict.edr}}}
\examples{
require(EDR)
\dontrun{demo(edr_ex4)}
}
\keyword{ smooth  }% at least one, from doc/KEYWORDS
\keyword{ regression }% __ONLY ONE__ keyword per line
