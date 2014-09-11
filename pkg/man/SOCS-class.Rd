\name{SOCS-class}
\Rdversion{1.1}
\docType{class}
\alias{SOCS-class}

\title{Class \code{"SOCS"}}

\description{
  Nesterov-Todd scalings and scaled Lagrange-multipliers for SOC
  constraints. 
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("SOCS",
  ...)} or by using \code{nnos()}.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{Object of class \code{"list"} with elements
      \code{beta} and \code{v} for the computation of the Nesterov-Todd
      scaling and its inverse and the scaled variable \code{lambda}.}
  }
}

\examples{
showClass("SOCS")
}
