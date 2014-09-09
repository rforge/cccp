\name{PSDC-class}
\Rdversion{1.1}
\docType{class}
\alias{PSDC-class}
\alias{psdc}

\title{Class \code{"PSDC"}}

\description{
  Definition for constraints with respect to the positive semidefinite cone.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PSDC",
  ...)}, or by invoking the constructor function \code{psdc(Flist, F0)}.
}

\section{Slots}{
  \describe{
    \item{\code{G}:}{Object of class \code{"Matrix"}: A \eqn{(p(p + 1) /
      2 \times n)} matrix containing the coefficients pertinent to the
      variable \eqn{x}.}
    \item{\code{h}:}{Object of class \code{"PSDV"}: Object containing
      the coefficients pertinent to the righthand-side term.}
    \item{\code{dims}:}{Object of class \code{"integer"}: The dimension
      of the cone.}
    \item{\code{vclass}:}{Object of class \code{"character"}: The class
      name of variables associated with the second-order cone.}
  }
}

\details{
  A psd-cone constraint is given as \eqn{\sum_{i = 1}^n x_i F_i \le
  F_0}. The matrix \eqn{G} is created as \eqn{G = [\textrm{vech}(F_1) |
  \ldots | \textrm{vech}(F_n)]} and the vector \eqn{h} is constructed as
  \eqn{h = [\textrm{vech}(F_0)]}.  
}

\examples{
showClass("PSDC")
}

\keyword{classes}
