\name{umsa}
\alias{umsa}
\alias{umsa,LNLV-method}
\alias{umsa,SOCV-method}
\alias{umsa,PSDV-method}

\title{
  Adjusting cone variable by maximum step size
}

\description{
  Returns for the primal and/or dual slack variables the maximum step
  size for adjusting the initial solution with respect to the cone
  constraint(s). 
}

\usage{
\S4method{umsa}{LNLV}(u, alpha, init = TRUE)  
\S4method{umsa}{SOCV}(u, alpha, init = TRUE)  
\S4method{umsa}{PSDV}(u, alpha, init = TRUE, sigma = NULL, lambda = NULL)  
}

\arguments{
  \item{u}{An object of either S4-class \code{LNLV}, \code{SOCV} or
    \code{PSDV}.}
  \item{alpha}{\code{numeric}, the maximim step size, ordinarily returned by a
    call to \code{maxs()}.}  
  \item{init}{\code{logical}, whether step size should be applied to
    initial solution.}
  \item{sigma}{\code{vector}, the eigen values of \code{u}.}
  \item{lambda}{An object of S4-class \code{PSDV} of the
    Lagrange-multipliers.}
}

\details{
  The initial values of the primal slack variable \eqn{\hat{s}} is given
  as \eqn{\hat{s} = -z} for \eqn{\alpha_p < 0} and \eqn{-z + (1 +
  \alpha_p) e} otherwise, whereby \eqn{\alpha_p = inf(\alpha | -z +
  \alpha e >= 0)}.\cr 
  The initial values of the dual slack variable \eqn{\hat{z}} is given
  as \eqn{\hat{z} = z} for \eqn{\alpha_d < 0} and \eqn{z + (1 +
  \alpha_d) e} otherwise, whereby \eqn{\alpha_d = inf(\alpha | z +
  \alpha e >= 0)}.\cr 
}

\value{
  An (adjusted) object of the same S4-class as object \code{u}.
}

\keyword{optimize}
