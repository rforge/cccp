\name{uvec}
\alias{uvec}
\alias{uvec,LNLV-method}
\alias{uvec,SOCV-method}
\alias{uvec,PSDV-method}

\title{
  Stacking of a symmetric matrix
}

\description{
  This function stacks the lower triangular and the main-diagonal
  elements of a matrix to a (scaled) vector. Hereby, the
  off-diagonal matrix elements are multiplied by \eqn{\sqrt{2}}.  
}

\usage{
\S4method{uvec}{LNLV}(u)  
\S4method{uvec}{SOCV}(u)  
\S4method{uvec}{PSDV}(u)  
}

\arguments{
  \item{u}{An object belonging to a class in the usage-section.}
}

\details{
  The scaling of the off-diagonal matrix elements ensures that inner
  products are preserved: \eqn{\textrm{tr}(UV) =
  \textrm{vech}(U)'\textrm{vech}(V)} for all U, V.  
}

\value{
  An object of the same class in packed storage.
}

\keyword{Optimize}
