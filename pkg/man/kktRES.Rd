\name{kktRES}
\alias{kktRES}
\alias{kktRES,DEFNL,PDV,PDV-method}
\alias{kktRES,DEFLP,PDV,PDV-method}
\alias{kktRES,DEFQP,PDV,PDV-method}

\title{
  Generic for computing the residuals of a KKT-system
}

\description{
  This function returns the residuals of a KKT-system.
}

\usage{
\S4method{kktRES}{DEFNL,PDV,PDV}(cpd, u, v, W, sdv)  
\S4method{kktRES}{DEFLP,PDV,PDV}(cpd, u, v, W, sdv)  
\S4method{kktRES}{DEFQP,PDV,PDV}(cpd, u, v, W)  
}

\arguments{
\item{cpd}{An object belonging to the S4-class union \code{CPD}.}
\item{u}{An object of S4-class \code{PDV}} 
\item{v}{An object of S4-class \code{PDV}} 
\item{W}{A \code{list} object with elements of S4-class \code{NTS}.}
\item{sdv}{A named vector with elements \sQuote{\code{dg}},
  \sQuote{\code{dgi}}, \sQuote{\code{lambdag}}: The values related to the
  self-dual variables \eqn{\tau}{tau} and \eqn{\kappa}{kappa}.}
}

\value{
  An object of S4-class \code{PDV}.
}

\keyword{optimize}
