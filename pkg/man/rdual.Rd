\name{rdual}
\alias{rdual}
\alias{rdual,PDV,DEFLP-method}
\alias{rdual,PDV,DEFQP-method}

\title{
  Residuals of dual constraints
}

\description{
  Returns the vector of dual residuals evaluated at the point \code{x,
    y, s, z}.
  
}

\usage{
\S4method{rdual}{PDV,DEFLP}(pdv, cpd)  
\S4method{rdual}{PDV,DEFQP}(pdv, cpd)  
}

\arguments{
  \item{pdv}{An object of S4-class \code{PDV} for computing the dual residuals.}
  \item{cpd}{An object belonging to the S4-class union \code{CPD}.}
  
}

\details{
  For problems with equality/inequality constraints, the dual
  residuals are given by:
  \deqn{r_\textrm{dual} = \nabla f(x) + A'y + G'z}
}

\value{
  \code{vector}, the residuals of the dual constraints evaluated at
  \code{x, y, s, z}.
}

\seealso{
  \code{\linkS4class{DEFLP}}, \code{\linkS4class{DEFQP}}, \code{\linkS4class{PDV}}
}

\keyword{optimize}
