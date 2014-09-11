\name{NLFV-class}
\Rdversion{1.1}
\docType{class}
\alias{NLFV-class}

\title{Class \code{"NLFV"}}

\description{
  Variable with respect to the nonlinear functions.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NLFV",
  ...)} or by using \code{nlfv(u)}.
}

\section{Slots}{
  \describe{
    \item{\code{u}:}{Object of class \code{"Matrix"}: A vector
      representing a point with respect to the nonnegative orthant cone.}
    \item{\code{dims}:}{Object of class \code{"integer"}: The dimension
      of the cone, \emph{i.e.}, the length of \code{u}.}
  }
}

\examples{
showClass("NLFV")
}
