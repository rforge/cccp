\name{SOCV-class}
\Rdversion{1.1}
\docType{class}
\alias{SOCV-class}

\title{Class \code{"SOCV"}}

\description{
  Variable with respect to the second-order cone.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SOCV",
  ...)}, or by using \code{socv(u)}.
}

\section{Slots}{
  \describe{
    \item{\code{u}:}{Object of class \code{"matrix"}: A vector
      for a variable with respect to the second-order cone.} 
    \item{\code{dims}:}{Object of class \code{"integer"}: The dimension
      of the cone.}
  }
}

\examples{
showClass("SOCV")
}

\keyword{classes}
