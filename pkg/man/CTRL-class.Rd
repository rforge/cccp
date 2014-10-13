\name{CTRL-class}
\Rdversion{1.1}
\docType{class}
\alias{CTRL-class}
\alias{validCTRL}

\title{Class \code{"CTRL"}}

\description{
  Class definition for control parameters utilized in optimization routines.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CTRL", ...)},
or more conveniently by invoking \code{ctrl()}.
}

\section{Slots}{
  \describe{
    \item{\code{maxiters}:}{Object of class \code{"integer"}: The
      maximum count of iterations.}
    \item{\code{abstol}:}{Object of class \code{"numeric"}: The absolute
      convergence level.}
    \item{\code{reltol}:}{Object of class \code{"numeric"}: The relative
      convergence level.}
    \item{\code{feastol}:}{Object of class \code{"numeric"}: The
      feasible convergence level.}
    \item{\code{refine}:}{Object of class \code{"logical"}: Whether
      iterative refinement should be applied.}
    \item{\code{stepadj}:}{Object of class \code{"numeric"}: Step size
      adjustment in combined step.}
    \item{\code{alpha}:}{Object of class \code{"numeric"}: Parameter in
      backtracking line search.}
    \item{\code{beta}:}{Object of class \code{"numeric"}: Parameter in
      backtracking line search.}
    \item{\code{maxreliter}:}{Object of class \code{"integer"}:
      Parameter in backtracking line search.}
    \item{\code{trace}:}{Object of class \code{"logical"}: \code{TRUE}
      if solver's progress should be returned during iterations,
      otherwise \code{FALSE}.}
    \item{\code{method}:}{Object of class \code{"character"}: The method
      for solving the KKT-system.}
  }
}

\details{
  The values for the slots are/can be validated by the function \code{validCTRL}.
}

\keyword{classes}
