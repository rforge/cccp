\name{PSDV-class}
\Rdversion{1.1}
\docType{class}
\alias{PSDV-class}
\alias{psdv}

\title{Class \code{"PSDV"}}

\description{
  Variable with respect to positive semidefinite cone.
}

\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("PSDV",
    ...)}, or by calling \code{psdv(U, degree = ncol(U)}.
}

\section{Slots}{
  \describe{
    \item{\code{u}:}{Object of class \code{"Matrix"}: A matrix
      representing a point with respect to the positive semidefinte cone
      in \code{vech}-notation.}
    \item{\code{dims}:}{Object of class \code{"integer"}: The degree
      of the cone.}
  }
}

\examples{
showClass("PSDV")
}
