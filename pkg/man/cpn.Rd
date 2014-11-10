\name{cpn}
\alias{cpn}

\title{
  Convex program with nonlinear constraints
}

\description{
  This function minimizes a linear objective with nonlinear and cone
  constraints. 
}

\usage{
cpn(x0, f0, nlfList = list(), cList = list(), A = NULL, b = NULL,
    optctrl = ctrl())
}

\arguments{
  \item{x0}{An object of class \code{vector} with dimension \eqn{n
      \times 1} for the initial values. The point \code{x0} must be in
    the domain of the nonlinear constraints.}
  \item{f0}{\code{function}: the convex and twice-differentiable
    objective function (its first argument must be \sQuote{\code{x}}).}
  \item{nlfList}{A \code{list} object containing the nonlinear
    constraints as its elements. The functions have to be specified with
    \code{x} as their first argument and must be casted in implicit
    form, \emph{i.e.} \eqn{f(x) \le 0}.}
  \item{cList}{A \code{list} object containing the cone
    constraints. Elements must be of either S4-class \code{NNOC}, or
    \code{SOCC}, or \code{PSDC} or an empty list in case of no
    inequality constraints.}
  \item{A}{An object of class \code{matrix} with dimension \eqn{p \times
      n} or \code{NULL} for problems without equality constraints.}
  \item{b}{An object of class \code{vector} with dimension \eqn{p \times
      1} or \code{NULL} for problems without equality constraints.}
  \item{optctrl}{An object of S4-class \code{CTRL}.}
}

\value{
  Returns an object of S4-class \code{CPS}.
}

\keyword{optimize}
