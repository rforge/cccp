\name{uslb}
\alias{uslb}
\alias{uslb,LNLV,LNLV-method}
\alias{uslb,SOCV,SOCV-method}
\alias{uslb,PSDV,PSDV-method}

\title{
  Log-Barrier scaling of cone variables.
}

\description{
  This function applies the scaling of a cone variables by means of the
  Hessian of the Log-Barrier function. 
}

\usage{
\S4method{uslb}{LNLV,LNLV}(u, lambda, inv = FALSE)  
\S4method{uslb}{SOCV,SOCV}(u, lambda, inv = FALSE)  
\S4method{uslb}{PSDV,PSDV}(u, lambda, inv = FALSE)  
}

\arguments{
\item{u}{An object of either S4-class \code{LNLV} or \code{SOCV} or \code{PSDV}.}
\item{lambda}{An object of the same S4-class as \code{u}.}
\item{inv}{\code{logical}, whether the inverse of the scaling matrix
  should be used.}
}

\value{
  An S4-object of the corresponding class for cone variables,
  \emph{e.g.}, \code{LNLV}.
}

\keyword{optimize}
