\name{uprd}
\alias{uprd}
\alias{uprd,LNLV,LNLV-method}
\alias{uprd,LNLV,missing-method}
\alias{uprd,SOCV,SOCV-method}
\alias{uprd,SOCV,missing-method}
\alias{uprd,PSDV,PSDV-method}
\alias{uprd,PSDV,missing-method}

\title{
  Product of cone-variable(s)
}

\description{
  This function returns the product of two points with respect to
  either the nonnegative orthant, the second-order or the semidefinite
  cone.    
}

\usage{
\S4method{uprd}{LNLV,LNLV}(u, v)  
\S4method{uprd}{LNLV,missing}(u, v)  
\S4method{uprd}{SOCV,SOCV}(u, v)  
\S4method{uprd}{SOCV,missing}(u, v)  
\S4method{uprd}{PSDV,PSDV}(u, v)  
\S4method{uprd}{PSDV,missing}(u, v)  
}

\arguments{
  \item{u}{A vector or an object of either S4-class \code{LNLV}, or
  \code{SOCV}, or \code{PSDV}.} 
  \item{v}{A vector or an object of either S4-class \code{LNLV}, or
  \code{SOCV}, or \code{PSDV}.} 
}

\details{
  For the nonnegative orthant cone, the dot-product is defined as:
  \eqn{u \circ v = (u_1 * v_1, \ldots, u_p * v_p)} for \eqn{C_k =
  R^p_+}.\cr 
  For the second-order cone, the dot-product is defined as:
  \eqn{u \circ v = (u'v, u_0 v_1  + v_0 u_1)}.\cr 
  For the positive semidefinite cone, the dot-product is defined as:
  \eqn{u \circ v = 0.5 *
  \textrm{vech}(\textrm{unvech}(u)\textrm{unvech}(v) +
  \textrm{unvech}(v)\textrm{unvech}(u))}.\cr
}

\value{
  An S4-object of the corresponding class for cone variables,
  \emph{e.g.}, \code{LNLV}.
}

\keyword{optimize}
