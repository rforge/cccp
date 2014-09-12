\name{cps}
\alias{cps}
\alias{cps,DEFNL-method}
\alias{cps,DEFLP-method}
\alias{cps,DEFQP-method}

\title{
  Solving a convex program
}

\description{
  This function returns an optimal point for a cone constraint convex
  program. The cpd is a member of the class union \code{CPD}. 
}

\usage{
\S4method{cps}{DEFNL}(cpd)  
\S4method{cps}{DEFLP}(cpd)  
\S4method{cps}{DEFQP}(cpd)  
}

\arguments{
  \item{cpd}{An object belonging to the class union \code{CPD}.}
}

\value{
  An object of S4-class \code{CPS}. 
}

\keyword{optimize}
