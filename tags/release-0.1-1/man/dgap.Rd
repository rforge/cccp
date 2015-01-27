\name{dgap}
\alias{dgap}
\alias{dgap,PDV-method}

\title{
  Generic for computing the duality gap
}

\description{
  Computes the duality gap defined as the dot product between the slack
  variables in the primal and dual formulation: \eqn{s' z}.
}

\usage{
\S4method{dgap}{PDV}(pdv)  
}

\arguments{
  \item{pdv}{An object of S4-class \code{PDV}}
}

\value{
  \code{numeric}
}

\keyword{optimize}
