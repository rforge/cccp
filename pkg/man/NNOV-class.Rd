\name{NNOV-class}
\Rdversion{1.1}
\docType{class}
\alias{NNOV-class}

\title{Class \code{"NNOV"}}

\description{
  Variable with respect to the non-negative orthant cone.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NNOV",
  ...)} or by using \code{nnov(u)}.
}

\section{Slots}{
  \describe{
    \item{\code{u}:}{Object of class \code{"matrix"}: A vector
      representing a point with respect to the nonnegative orthant cone.}
    \item{\code{dims}:}{Object of class \code{"integer"}: The dimension
      of the cone, \emph{i.e.}, the length of \code{u}.}
  }
}

\examples{
showClass("NNOV")
}
