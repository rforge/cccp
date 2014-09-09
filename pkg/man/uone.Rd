\name{uone}
\alias{uone}
\alias{uone,NNOV-method}
\alias{uone,SOCV-method}
\alias{uone,PSDV-method}

\title{
  One-element of a cone-variable
}

\description{
  This function returns the one-element of a variable with respect
  to either the nonnegative orthant, the second-order or the
  semidefinite cone. 
}

\usage{
\S4method{uone}{NNOV}(u)  
\S4method{uone}{SOCV}(u)  
\S4method{uone}{PSDV}(u)  
}

\arguments{
  \item{u}{A vector or an object of either S4-class \code{NNOV}, or
  \code{SOCV}, or \code{PSDV}.} 
}

\details{
  For the nonnegative orthant cone, the one-element \eqn{\bm{e}} is
  defined as the one-vector, \emph{i.e.} \eqn{e = (1, 1, \ldots, 1)}.\cr 
  For the second-order cone, the uone-element is definded as \eqn{e =
    (1, 0, \ldots, 0)}.\cr
  For the psoitive semidefinite cone, the uone-element is definded as
  \eqn{e = \textrm{vech}(I_p)}.
}

\value{
  An S4-object of the corresponding class for cone variables,
  \emph{e.g.}, \code{NNOV}.
}

\keyword{optimize}
