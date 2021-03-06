\name{rprim}
\alias{rprim}
\alias{rprim,PDV,CPD-method}

\title{
  Residuals of primal constraints
}

\description{
  Returns the vector of primal residuals evaluated at the point \code{x,
    y, s, z}.
}

\usage{
\S4method{rprim}{PDV,CPD}(pdv, cpd)  
}

\arguments{
  \item{pdv}{An object of S4-class \code{PDV} for computing the dual residuals.}
  \item{cpd}{An object belonging to the S4-class union \code{CPD}.}
  
}

\details{
  For problems with equality constraints, the primal residuals are given by:
  \deqn{r_\textrm{prim} = b - A'x}.
}

\value{
  \code{vector}, the residuals of the dual constraints evaluated at
  \code{x, y, s, z}.
}

\seealso{
  \code{\linkS4class{DEFLP}}, \code{\linkS4class{DEFQP}},
  \code{\linkS4class{DEFNL}}, \code{\linkS4class{DEFCP}},
  \code{\linkS4class{PDV}}
}

\keyword{optimize}
