\name{umat}
\alias{umat}
\alias{umat,NNOV-method}
\alias{umat,SOCV-method}
\alias{umat,PSDV-method}

\title{
  Inverse vech-operation
}

\description{
  This method converts a vector (one-column matrix) of an
  S4-class object \code{PSDV} in its (symmetric)
  matrix-representation. The off-diagonal entries of the resulting
  Matrix are divided by \eqn{\sqrt{2}}.
}

\usage{
umat(u)
\S4method{umat}{NNOV}(u)  
\S4method{umat}{SOCV}(u)  
\S4method{umat}{PSDV}(u)  
}

\arguments{
  \item{u}{An object belonging to a class in the usage-section.}
}

\details{
  If \code{u} is a one-column matrix with row dimension equal to \eqn{p
  * (p + 1) / 2}, then the scaling of the off-diagonal elements ensures
  that \eqn{u'v = \textrm{tr}(\textrm{umat}(u) \textrm{umat}(v))}.   
}

\value{
  An object of the same S4-class as its argument.
}

\keyword{optimize}
