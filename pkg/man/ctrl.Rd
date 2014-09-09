\name{ctrl}
\alias{ctrl}

\title{
  Creating objects of S4-class CTRL
}

\description{
  This function creates an object of S4-class \code{CTRL} which
  contains optimization parameters, \emph{e.g.} the maximum number of
  iterations. 
}

\usage{
ctrl(maxiters = 100L, abstol = 1e-07, reltol = 1e-06,
     feastol = 1e-07, refine = FALSE, stepadj = 0.95, trace = TRUE, 
     method = c("solve"))
}

\arguments{
  \item{maxiters}{\code{integer}, the maximum count of iterations.}
  \item{abstol}{\code{numeric}, the absolute level for convergence to be
    achieved.}
  \item{reltol}{\code{numeric}, the relative level for convergence to be
    achieved.}
  \item{feastol}{\code{numeric}, the feasable level for convergence to be
    achieved.}
  \item{refine}{\code{logical}, if \code{FALSE} (the default), no
    iterative refinement is applied.}
  \item{stepadj}{\code{numeric}, step size adjustment in combined step.}
  \item{trace}{\code{logical}, if \code{TRUE} (the default), the
    solver's progress during the iterations is shown.}
  \item{method}{\code{"character"}, the method for solving the KKT-system.}
}

\value{
  An object of S4-class \code{CTRL}.
}

\note{
  Either \code{abstol} or \code{reltol} can be set to a negative real
  number. \code{feastol} must be greater than zero.
}

\seealso{
  \code{\linkS4class{CTRL}}
}

\keyword{optimize}
