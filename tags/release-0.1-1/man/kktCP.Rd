\name{kktCP}
\alias{kktCP}

\title{
  Function for computing a solution of a KKT-system for convex programs
}

\description{
  This function returns the solution of a KKT-system for convex
  programs, \emph{i.e.} S4-objects of class \code{DEFCP}.
}

\usage{
kktCP(cpd)  
}

\arguments{
\item{cpd}{An object belonging to the S4-class \code{DEFCP}.}
}

\value{
  An object of S4-class \code{KKTSLV}.
}

\keyword{optimize}
