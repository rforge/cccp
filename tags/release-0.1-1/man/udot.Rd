\name{udot}
\alias{udot}
\alias{udot,numeric,numeric-method}
\alias{udot,numeric,missing-method}
\alias{udot,LNLV,LNLV-method}
\alias{udot,LNLV,missing-method}
\alias{udot,SOCV,SOCV-method}
\alias{udot,SOCV,missing-method}
\alias{udot,PSDV,PSDV-method}
\alias{udot,PSDV,missing-method}

\title{
  Inner-product between cone-variables
}

\description{
  This function returns the inner-product between two points with
  respect to either the nonnegative orthant, the second-order or the
  semidefinite cone.    
}

\usage{
\S4method{udot}{numeric,numeric}(u, v)  
\S4method{udot}{numeric,missing}(u, v)  
\S4method{udot}{LNLV,LNLV}(u, v)  
\S4method{udot}{LNLV,missing}(u, v)  
\S4method{udot}{SOCV,SOCV}(u, v)  
\S4method{udot}{SOCV,missing}(u, v)  
\S4method{udot}{PSDV,PSDV}(u, v)  
\S4method{udot}{PSDV,missing}(u, v)  
}

\arguments{
  \item{u}{A vector or an object of either S4-class \code{LNLV}, or
  \code{SOCV}, or \code{PSDV}.} 
  \item{v}{An object of the same class as \code{u}.} 
}


\value{
  \code{numeric}.
}

\keyword{optimize}
