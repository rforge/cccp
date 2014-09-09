\name{ntsu}
\alias{ntsu}
\alias{ntsu,NNOS,NNOV,NNOV-method}
\alias{ntsu,SOCS,SOCV,SOCV-method}
\alias{ntsu,PSDS,PSDV,PSDV-method}

\title{
  Updating Nesterov-Todd scalings and Lagrange-multipliers
}

\description{
  This method updates the Nesterov-Todd scaling matrix, its inverse and
  the Lagrange multipliers of points for the primal and dual
  slack-variables with respect to a self-dual cone.
}

\usage{
\S4method{ntsu}{NNOS,NNOV,NNOV}(W, s, z)  
\S4method{ntsu}{SOCS,SOCV,SOCV}(W, s, z)  
\S4method{ntsu}{PSDS,PSDV,PSDV}(W, s, z)  
}

\arguments{
  \item{W}{An object of either S4-class \code{NNOS} or \code{SOCS} or \code{PSDS}.} 
  \item{s}{An object of either S4-class \code{NNOV} or \code{SOCV} or \code{PSDV}.} 
  \item{z}{An object of either S4-class \code{NNOV} or \code{SOCV} or \code{PSDV}.}
}

\seealso{
  \code{\linkS4class{NNOV}}, \code{\linkS4class{NNOS}}, 
  \code{\linkS4class{SOCV}}, \code{\linkS4class{SOCS}},
  \code{\linkS4class{PSDV}}, \code{\linkS4class{PSDS}} 
}

\value{
  An object of either S4-class \code{NNOS} or \code{SOCS} or \code{PSDS}.
}

\keyword{Optimize}
