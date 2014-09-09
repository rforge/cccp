\name{jdot}
\alias{jdot}
\alias{jdot,numeric,numeric-method}
\alias{jdot,SOCV,SOCV-method}
\alias{jdot,SOCV,missing-method}

\title{
  Scaling constant in hyperbolic Householder transformation
}

\description{
  This function returns the scalar \eqn{u' J v}, whereby the matrix
  \eqn{J} is defined as \eqn{J = [1, 0; 0, -I]}. 
}

\usage{
\S4method{jdot}{numeric,numeric}(u, v)  
\S4method{jdot}{SOCV,SOCV}(u, v)  
\S4method{jdot}{SOCV,missing}(u, v)  
}

\arguments{
  \item{u}{An object of either S4-class \code{SOCV} or \code{numeric}.} 
  \item{v}{An object of S4-class \code{SOCV}.} 
}


\value{
  \code{numeric}.
}

\keyword{optimize}
