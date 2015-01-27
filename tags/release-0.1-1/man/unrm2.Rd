\name{unrm2}
\alias{unrm2}
\alias{unrm2,numeric-method}
\alias{unrm2,LNLV-method}
\alias{unrm2,SOCV-method}
\alias{unrm2,PSDV-method}

\title{
  Norm of a cone-variable
}

\description{
  This function returns the norm of a variable with respect
  to a cone variable. 
}

\usage{
\S4method{unrm2}{numeric}(u)  
\S4method{unrm2}{LNLV}(u)  
\S4method{unrm2}{SOCV}(u)  
\S4method{unrm2}{PSDV}(u)  
}

\arguments{
  \item{u}{A vector or an object of either S4-class \code{LNLV},
  \code{SOCV}, or \code{PSDV}.} 
}

\value{
  \code{numeric}.
}

\keyword{optimize}
