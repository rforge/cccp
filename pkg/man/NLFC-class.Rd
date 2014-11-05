\name{NLFC-class}
\Rdversion{1.1}
\docType{class}
\alias{NLFC-class}
\alias{nlfc}

\title{Class \code{"NLFC"}}

\description{
  Definition for nonlinear constraints.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NLFC",
  ...)}, or by invoking the constructor function \code{nlfc(G, h)}.
}

\section{Slots}{
  \describe{
    \item{\code{G}:}{Object of class \code{"matrix"}: A \eqn{(m \times n)} matrix
      containing the coefficients of the lefthand-side linear inequality
      constraints.}
    \item{\code{h}:}{Object of class \code{"NLFV"}: A \eqn{(m \times 1)} vector
      containing the coefficients of the righthand-side linear inequality
      constraints as slot \code{u}.}
    \item{\code{dims}:}{Object of class \code{"integer"}: The dimension
      of the cone.}
    \item{\code{vclass}:}{Object of class \code{"character"}: The class
      name of variables associated with the NLF-cone.}
  }
}

\examples{
showClass("NLFC")
}

\keyword{classes}
