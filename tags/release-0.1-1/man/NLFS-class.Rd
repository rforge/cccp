\name{NLFS-class}
\Rdversion{1.1}
\docType{class}
\alias{NLFS-class}

\title{Class \code{"NLFS"}}

\description{
  Nesterov-Todd scalings and scaled Lagrange-multipliers for nonlinear
  constraints. 
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("NLFS",
  ...)} or by using \code{nnos()}.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{Object of class \code{"list"} with elements
      \code{d} and \code{di} for the Nesterov-Todd scaling and its
      inverse and the scaled variable \code{lambda}.}
  }
}

\examples{
showClass("NLFS")
}
