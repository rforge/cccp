\name{cpd}
\alias{cpd}

\title{
  Creating a member object of the class union \code{CPD}
}

\description{
  This function returns an object containing the definition of a convex
  program. The object is a member of the class union \code{CPD}. 
}

\usage{
cpd(P, q, A = NULL, b = NULL, cList = list(), optctrl = ctrl())
}

\arguments{
  \item{P}{An object of class \code{Matrix} with dimension \eqn{n \times n}.}
  \item{q}{An object of class \code{vector} with dimension \eqn{n \times 1}.}
  \item{A}{An object of class \code{Matrix} with dimension \eqn{p \times
      n} or \code{NULL} for problems without equality constraints.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times
      1} or \code{NULL} for problems without equality constraints.}
  \item{cList}{A \code{list} object containing the cone
    constraints. Elements must be of either S4-class \code{NNOC}, or
    \code{SOCC}, or \code{PSDC} or an empty list in case of no
    inequality constraints.}
  \item{optctrl}{An object of S4-class \code{CTRL}.}
}

\value{
  An object belonging to the class union \code{CPD}. 
}

\keyword{optimize}
