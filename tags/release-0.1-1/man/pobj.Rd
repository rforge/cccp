\name{pobj}
\alias{pobj}
\alias{pobj,numeric,DEFCP-method}
\alias{pobj,PDV,DEFCP-method}
\alias{pobj,numeric,DEFNL-method}
\alias{pobj,PDV,DEFNL-method}
\alias{pobj,numeric,DEFLP-method}
\alias{pobj,PDV,DEFLP-method}
\alias{pobj,numeric,DEFQP-method}
\alias{pobj,PDV,DEFQP-method}

\title{
  Value of the primal objective
}

\description{
  Returns the value of the primal objective function evaluated at the
  point \code{x}.  
}

\usage{
\S4method{pobj}{numeric,DEFCP}(pdv, cpd)  
\S4method{pobj}{PDV,DEFCP}(pdv, cpd)  
\S4method{pobj}{numeric,DEFNL}(pdv, cpd)  
\S4method{pobj}{PDV,DEFNL}(pdv, cpd)  
\S4method{pobj}{numeric,DEFLP}(pdv, cpd)  
\S4method{pobj}{PDV,DEFLP}(pdv, cpd)  
\S4method{pobj}{numeric,DEFQP}(pdv, cpd)  
\S4method{pobj}{PDV,DEFQP}(pdv, cpd)  
}

\arguments{
  \item{pdv}{An object of either S4-class \code{PDV} or \code{numeric}
    for evaluating the objective function.}
  \item{cpd}{An object belonging to the S4-class union \code{CPD}.}
}

\details{
  For linear problems the objective is defined as:
  \deqn{f(x) = q'x}
  
  For quadratic problems the objective is defined as:
  \deqn{f(x) = 1/ 2 x' P x + q'x}
}

\value{
  \code{numeric}, the primal objective evaluated at the point \eqn{x}.
}

\seealso{
    \code{\linkS4class{DEFLP}}, \code{\linkS4class{DEFQP}},
    \code{\linkS4class{DEFNL}}, \code{\linkS4class{DEFCP}},
    \code{\linkS4class{PDV}}
}

\keyword{optimize}
