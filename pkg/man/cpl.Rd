\name{cpl}
\alias{cpl}

\title{
  Linear program with nonlinear constraints
}

\description{
  This function minimizes a linear objective with nonlinear and cone
  constraints. 
}

\usage{
cpl(x0, q, nlfList = list(), cList = list(), A = NULL, b = NULL,
    optctrl = ctrl())
}

\arguments{
  \item{x0}{An object of class \code{vector} with dimension \eqn{n
      \times 1} for the initial values. The point \code{x0} must be in
    the domain of the nonlinear constraints.}
  \item{q}{\code{vector} of length \eqn{n} for the coefficients in the
    objective.}
  \item{nlfList}{A \code{list} object containing the nonlinear
    constraints as its elements. The functions have to be specified with
    \code{x} as their first argument and must be casted in implicit
    form, \emph{i.e.} \eqn{f(x) \le 0}.}
  \item{cList}{A \code{list} object containing the cone
    constraints. Elements must be of either S4-class \code{NNOC}, or
    \code{SOCC}, or \code{PSDC} or an empty list in case of no
    inequality constraints.}
  \item{A}{An object of class \code{Matrix} with dimension \eqn{p \times
      n} or \code{NULL} for problems without equality constraints.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times
      1} or \code{NULL} for problems without equality constraints.}
  \item{optctrl}{An object of S4-class \code{CTRL}.}
}

\value{
  Currently, aAn object of S4-class \code{DEFNL}.
}

\keyword{optimize}
