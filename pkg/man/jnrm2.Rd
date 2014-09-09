\name{jnrm2}
\alias{jnrm2}
\alias{jnrm2,SOCV-method}

\title{
  Norm of a scaling-point for SOC-variables
}

\description{
  This function returns the square-root of \eqn{u ' J u}, whereby \eqn{J
  = [1, 0; 0, -I]}. 
}

\usage{
\S4method{jnrm2}{SOCV}(u)  
}

\arguments{
  \item{u}{An object of S4-class \code{SOCV}.} 
}

\value{
  \code{numeric}.
}

\keyword{optimize}
