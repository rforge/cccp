\name{kktSOL}
\alias{kktSOL}
\alias{kktSOL,DEFLP-method}
\alias{kktSOL,DEFQP-method}

\title{
  Generic for computing a solution of a KKT-system
}

\description{
  This function returns the solution of a KKT-system.
}

\usage{
\S4method{kktSOL}{DEFLP}(cpd, SolKkt2, SolKkt1, W, WhL, sdv, kktslv, refine = FALSE)  
\S4method{kktSOL}{DEFQP}(cpd, SolKkt, W, kktslv, refine = FALSE)  
}

\arguments{
\item{cpd}{An object belonging to the S4-class union \code{CPD}.}
\item{SolKkt}{An object of S4-class \code{PDV}: The variables for
  computing the RHS-values \code{x, y, z} used in the call to
  \code{kktslv\@f()}.} 
\item{SolKkt2}{An object of S4-class \code{PDV}: The variables for
  computing the RHS-values \code{x, y, z} used in the call to
  \code{kktslv\@f()}.} 
\item{SolKkt1}{An object of S4-class \code{PDV}: The solution of the
  first KKT-system.}
\item{W}{A \code{list} object with elements belonging to the S4-class
  union \code{NTS}.}
\item{WhL}{A \code{list} object with elements belonging to the
  S$-class union \code{CPV}: the scaled RHS-points of the inequality
  constraints.} 
\item{sdv}{A named vector with elements \sQuote{\code{dg}},
  \sQuote{\code{dgi}}, \sQuote{\code{lambdag}}: The values related to the
  self-dual variables \eqn{\tau}{tau} and \eqn{\kappa}{kappa}.}
\item{kktslv}{An object of S4-class \code{KKTSLV}.}
\item{refine}{\code{logical}, whether an iterativ refinement should be applied.}
}

\value{
  An object of S4-class \code{PDV}.
}

\keyword{optimize}
