\name{DEFLP-class}
\Rdversion{1.1}
\docType{class}
\alias{DEFLP-class}

\title{Class \code{"DEFLP"}}

\description{
  Definition of linear programs.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("DEFLP",
  ...)}, or by calling \code{cpd()}.
}

\section{Slots}{
  \describe{
    \item{\code{q}:}{Object of class \code{"vector"}: A \eqn{(n \times
	1)} vector of the coefficients in the objective function.}
    \item{\code{A}:}{Object of class \code{"Matrix"}: A \eqn{(p
	\times n)} matrix of the equality constraints (lhs).}
    \item{\code{b}:}{Object of class \code{"Matrix"}: A \eqn{(p
	\times 1)} matrix of the equality constraints (rhs).}
    \item{\code{cList}:}{Object of class \code{"list"}: The elements
      of \code{cList} are the cone constraints, each of either
      S4-class \code{NNOC}, or \code{SOCC}, or \code{PSDC}.}
    \item{\code{n}:}{Object of class \code{"integer"}: The count of
      variables in the objective function.}
    \item{\code{k}:}{Object of class \code{"integer"}: The count of
      cone constraints, \emph{i.e.}, the count of list elements in
      \code{cList}.}
    \item{\code{ctrl}:}{Object of class \code{"CTRL"}: The control
      parameters used in the IP-algorithm.}
    \item{\code{title}:}{Object of class \code{"character"}: A brief
      description of the quadratic program.}    
  }
}

\keyword{classes}
