\name{rcent}
\alias{rcent}
\alias{rcent,PDV,CPD-method}

\title{
  Residuals of centrality constraints
}

\description{
  Returns the vector of centrality residuals evaluated at the point \code{x, s}. 
}

\usage{
\S4method{rcent}{PDV,CPD}(pdv, cpd)  
}

\arguments{
  \item{pdv}{An object of S4-class \code{PDV} for computing the dual residuals.}
  \item{cpd}{An object of the S4-class union \code{CPD}.}
  
}


\details{
  For problems with inequality constraints, the cenrality residuals are given by:
  \deqn{r_\textrm{cent} = s + G x - h}
}

\value{
  \code{list}, the elements are the residuals of the centrality constraints evaluated at
  \code{x, s} for each cone constraint.
}

\seealso{
  \code{\linkS4class{DEFLP}}, \code{\linkS4class{DEFQP}}, \code{\linkS4class{PDV}}
}

\keyword{optimize}
