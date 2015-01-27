\name{ntsc}
\alias{ntsc}
\alias{ntsc,NLFV,NLFV-method}
\alias{ntsc,NLFC,missing-method}
\alias{ntsc,NNOV,NNOV-method}
\alias{ntsc,NNOC,missing-method}
\alias{ntsc,SOCV,SOCV-method}
\alias{ntsc,SOCC,missing-method}
\alias{ntsc,PSDV,PSDV-method}
\alias{ntsc,PSDC,missing-method}

\title{
  Nesterov-Todd scalings and scaled Lagrange-multipliers
}

\description{
  This method returns the Nesterov-Todd scaling matrix and its inverse
  of points for the primal and dual slack-variables with respect to a
  self-dual cone. An initial scaling is returned, if a cone constraint
  object is provided for the argument \code{s}.    
}

\usage{
\S4method{ntsc}{NLFV,NLFV}(s, z)  
\S4method{ntsc}{NLFC,missing}(s, z)  
\S4method{ntsc}{NNOV,NNOV}(s, z)  
\S4method{ntsc}{NNOC,missing}(s, z)  
\S4method{ntsc}{SOCV,SOCV}(s, z)  
\S4method{ntsc}{SOCC,missing}(s, z)  
\S4method{ntsc}{PSDV,PSDV}(s, z)  
\S4method{ntsc}{PSDC,missing}(s, z)  
}

\arguments{
  \item{s}{An object of S4-class \code{NNOV} or \code{NNOC}, or
  \code{SOCV} or \code{SOCC}, or \code{PSDV} or \code{PSDC}.} 
  \item{z}{An object of the same S4-class as \code{s}.}
}

\seealso{
  \code{\linkS4class{NLFV}}, \code{\linkS4class{NLFC}}, 
  \code{\linkS4class{NNOV}}, \code{\linkS4class{NNOC}}, 
  \code{\linkS4class{SOCV}}, \code{\linkS4class{SOCC}}, 
  \code{\linkS4class{PSDV}}, \code{\linkS4class{PSDC}}, 
}

\value{
  An object of either S4-class \code{NLFS}, or \code{NNOS}, or
  \code{SOCS}, or \code{PSDS}. 
}

\keyword{Optimize}
