\name{SOCC-class}
\Rdversion{1.1}
\docType{class}
\alias{SOCC-class}
\alias{socc}

\title{Class \code{"SOCC"}}

\description{
  Definition for constraints with respect to the second-order cone.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SOCC",
  ...)}, or by invoking the constructor function \code{socc(F, g, d, f)}.
}

\section{Slots}{
  \describe{
    \item{\code{G}:}{Object of class \code{"Matrix"}: A \eqn{(p \times n)} matrix
      containing the coefficients pertinent to the variable \eqn{x}.}
    \item{\code{h}:}{Object of class \code{"vector"}: A \eqn{(p \times 1)} vector
      containing the coefficients pertinent to the affine terms in the soc-constraints.}
    \item{\code{dims}:}{Object of class \code{"integer"}: The dimension
      of the cone.}
    \item{\code{vclass}:}{Object of class \code{"character"}: The class
      name of variables associated with the second-order cone.}
  }
}

\details{
  A second-order cone constraint is given as \eqn{|| Fx + g ||_2 \le d'x
  + f}. The matrix \eqn{G} is created as \eqn{G = [-d, -F]} and the
  vector \eqn{h} is constructed as \eqn{h = [f, g]}. 
}

\examples{
showClass("SOCC")
}

\keyword{classes}
