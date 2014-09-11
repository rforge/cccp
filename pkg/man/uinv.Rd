\name{uinv}
\alias{uinv}
\alias{uinv,LNLV,LNLV-method}
\alias{uinv,SOCV,SOCV-method}
\alias{uinv,PSDV,PSDV-method}

\title{
  Inverse of the product between two variables
}

\description{
  This function returns the inverse of the product between two
  points \eqn{u, v} with respect to either the nonnegative orthant, the
  second-order or the semidefinite cone. The following condition holds
  for any \eqn{v}: \eqn{u \dot (u \diamond v) = v}.  
}

\usage{
\S4method{uinv}{LNLV,LNLV}(u, v)  
\S4method{uinv}{SOCV,SOCV}(u, v)  
\S4method{uinv}{PSDV,PSDV}(u, v)  
}

\arguments{
  \item{u}{An object of either S4-class \code{LNLV}, or \code{SOCV}, or
  \code{PSDV}.}
  \item{v}{An object of the same S4-class as \code{u}.} 
}

\value{
  An S4-object of the same class as its inputs.
}

\keyword{optimize}
