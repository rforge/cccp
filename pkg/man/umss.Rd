\name{umss}
\alias{umss}
\alias{umss,NNOV-method}
\alias{umss,SOCV-method}
\alias{umss,PSDV-method}

\title{
  Computing the maximum step-size of cone variable(s)
}

\description{
  This function returns the maximum step size to be applied to the
  increments at the initial stage and after the affine and combined
  directions during the iterations.
}

\usage{
\S4method{umss}{NNOV}(u)  
\S4method{umss}{SOCV}(u)  
\S4method{umss}{PSDV}(u)  
}

\arguments{
  \item{u}{An object of either S4-class \code{NNOV}, \code{SOCV}, or
    \code{PSDV}.}
}

\value{
  A named \code{list} object with two elements: element \sQuote{\code{ms}} for the
  maximum step size and element \sQuote{\code{evd}} for the eigen value
  decomposition. For objects of either S4-class \code{NNOV} or
  \code{SOCV}, this element is set to \code{NULL}.
}

\keyword{optimize}
