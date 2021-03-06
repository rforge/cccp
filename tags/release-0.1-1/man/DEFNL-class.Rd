\name{DEFNL-class}
\Rdversion{1.1}
\docType{class}
\alias{DEFNL-class}

\title{Class \code{"DEFNL"}}

\description{
  Definition of linear programs with non-linear constraints.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DEFNL",
  ...)}, or by calling \code{cpd()}.
}

\section{Slots}{
  \describe{
    \item{\code{x0}:}{Object of class \code{"vector"}: A \eqn{(n \times
	1)} vector of initial values. This point must be in the domaine
      of the non-linear constraints.}
    \item{\code{q}:}{Object of class \code{"vector"}: A \eqn{(n \times
	1)} vector of the coefficients in the objective function.}
    \item{\code{nlfList}:}{Object of class \code{"list"}: A list with the
      scalar-valued non-linear constraint function(s) (first argument
      \code{x} as its elements.}
    \item{\code{nlgList}:}{Object of class \code{"list"}: A list with the
      gradient function(s) (first argument \code{x} as its elements.}
    \item{\code{nlhList}:}{Object of class \code{"list"}: A list with the
      Hessian function(s) (first argument \code{x} as its elements.}
    \item{\code{cList}:}{Object of class \code{"list"}: The elements are
      the constraints. The first element is pertinent to the nonlinear
      constraints (object of S4-class \code{NLFC}) and the remaining
      elements are the cone constraints, each of either S4-class
      \code{NNOC}, or \code{SOCC}, or \code{PSDC}, if applicable.}  
    \item{\code{A}:}{Object of class \code{"matrix"}: A \eqn{(p
	\times n)} matrix of the equality constraints (lhs).}
    \item{\code{b}:}{Object of class \code{"matrix"}: A \eqn{(p
	\times 1)} matrix of the equality constraints (rhs).}
    \item{\code{n}:}{Object of class \code{"integer"}: The count of
      variables in the objective function.}
    \item{\code{k}:}{Object of class \code{"integer"}: The count of
      cone constraints, \emph{i.e.}, the count of list elements in
      \code{cList}.}
    \item{\code{mnl}:}{Object of class \code{"integer"}: The count of
      non-linear constraints.}
    \item{\code{ctrl}:}{Object of class \code{"CTRL"}: The control
      parameters used in the IP-algorithm.}
    \item{\code{H}:}{Object of class \code{"matrix"}: A \eqn{(n \times
	n)} matrix containing the Hessian with respect to the nonlinear
      constraints.}
    \item{\code{title}:}{Object of class \code{"character"}: A brief
      description of the quadratic program.}    
  }
}

\keyword{classes}
