\name{PDV-class}
\Rdversion{1.1}
\docType{class}
\alias{PDV-class}

  
\title{Class \code{"PDV"}}

\description{
  Values of primal and dual variables in convex optimization problems.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PDV",
  ...)}.
}

\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"}: A \eqn{(N \times 1)}
      vector of the primal variable.}
    \item{\code{y}:}{Object of class \code{"numeric"}: A \eqn{(P \times 1)} vector
      of the dual variables associated with the equality constraints
      (dual problem).}
    \item{\code{s}:}{Object of class \code{"list"}: A list of
      length \eqn{K} with elements of the slack variables in the
      generalized inequality constraints in the primal problem. These
      are of either S4-class \code{NNOV}, \code{SOCV}, or \code{PSDV}.}
    \item{\code{z}:}{Object of class \code{"list"}: A list of
      length \eqn{K} with elements of the slack variables in the
      generalized inequality constraints in the dual problem. These
      are of either S4-class \code{NNOV}, \code{SOCV}, or \code{PSDV}.}    
    \item{\code{kappa}:}{Object of class \code{"numeric"}: A scalar
      related to solving the self-dual embedding of cone LPs.}
    \item{\code{tau}:}{Object of class \code{"numeric"}: A scalar
      related to solving the self-dual embedding of cone LPs.}
  }
}

\keyword{classes}
