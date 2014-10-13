\name{dobj}
\alias{dobj}
\alias{dobj,PDV,DEFCP-method}
\alias{dobj,PDV,DEFNL-method}
\alias{dobj,PDV,DEFLP-method}
\alias{dobj,PDV,DEFQP-method}

\title{
  Value of dual objective
}

\description{
  Returns the value of the dual objective function evaluated at the
  point \code{x}.  
}

\usage{
\S4method{dobj}{PDV,DEFCP}(pdv, cpd)  
\S4method{dobj}{PDV,DEFNL}(pdv, cpd)  
\S4method{dobj}{PDV,DEFLP}(pdv, cpd)  
\S4method{dobj}{PDV,DEFQP}(pdv, cpd)  
}

\arguments{
  \item{pdv}{An object of S4-class \code{PDV} for evaluating the objective function.}  
  \item{cpd}{An object belonging to the S4-class union \code{CPD}.}
}

\details{
  For linear programs with nonlinear inequality constraints the dual
  objective is defined as: 
  \deqn{L(x, y, znl, zl) = q'x + znl'* f(x) + zl'*(G x - h) + y'(A x - b)}
  
  For linear programs the dual objective is defined as:
  \deqn{L(y, z) = -h'z - b'y}

  For quadratic programs the dual objective is defined as:
  \deqn{L(x, y, z) = 1/ 2 x' P x + q'x + r + z'*(G x - h) + y'(A x - b)}
}

\value{
  \code{numeric}, the dual objective evaluated at the point \eqn{(x, y, z}.
}

\seealso{
  \code{\linkS4class{DEFLP}}, \code{\linkS4class{DEFQP}},
  \code{\linkS4class{DEFNL}}, \code{\linkS4class{PDV}}
}

\keyword{optimize}
