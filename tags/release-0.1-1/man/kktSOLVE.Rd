\name{kktSOLVE}
\alias{kktSOLVE}

\title{
  Solving KKT-system by LU decomposition
}

\description{
  This function returns a function with arguments \code{W} and
  \code{cpd} for solving the KKT-system by means of a LU decomposition;
  hereby the function \code{lu} from the package \pkg{Matrix} is employed.
}

\usage{
kktSOLVE(cpd)
}

\arguments{
  \item{cpd}{An object belonging to the S4-class union \code{CPD}.}
}

\value{
  A list with elements \code{f}: a function with arguments \code{W} and
  \code{cpd} and element \code{fitems}: a list with objects pertinent to
  \code{f}. \code{W} is a list of the scaling matrices and \code{cpd}
  belongs to the class union \code{CPD}. 
}

\seealso{
  \code{\linkS4class{CPD}}, \code{\linkS4class{NNOS}} 
}

\keyword{optimize}
