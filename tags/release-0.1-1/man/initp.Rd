\name{initp}
\alias{initp}
\alias{initp,NLFC-method}
\alias{initp,NNOC-method}
\alias{initp,SOCC-method}
\alias{initp,PSDC-method}

\title{
  Creating initial point for primal/dual slack variables
}

\description{
  These methods create an initial point for the slack variables.
}

\usage{
\S4method{initp}{NLFC}(object)  
\S4method{initp}{NNOC}(object)  
\S4method{initp}{SOCC}(object)  
\S4method{initp}{PSDC}(object)  
}

\arguments{
  \item{object}{An object of either S4-class \code{NLFC}, or
  \code{NNOC}, or \code{SOCC}, or \code{PSDC}.}
}

\value{
  An S4-object of the corresponding variable class, \emph{e.g.},
  \code{NLFV}.   
}

\keyword{optimize}
