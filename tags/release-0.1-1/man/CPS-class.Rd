\name{CPS-class}
\Rdversion{1.1}
\docType{class}
\alias{CPS-class}
\alias{show,CPS-method}

\title{Class \code{"CPS"}}

\description{
  Class for the solution of convex program.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("CPS", ...)} or by
invoking \code{cps()}.
}

\section{Slots}{
  \describe{
    \item{\code{x}:}{Object of class \code{"numeric"}: The values for the
      primal variables.}
    \item{\code{y}:}{Object of class \code{"numeric"}: The values
      of the dual variables associated with the equality constraints
      (primal formulation).}
    \item{\code{s}:}{Object of class \code{"list"}: The values f
      the dual variables associated with inequality constraints (primal
      formulation).}
    \item{\code{z}:}{Object of class \code{"list"}: The values
      of the dual variables associated with inequality constraints (dual
      formulation).} 
    \item{\code{pobj}:}{Object of class \code{"numeric"}: The value of
      the primal objective.}
    \item{\code{dobj}:}{Object of class \code{"numeric"}: The value
      of the dual objective.}
    \item{\code{dgap}:}{Object of class \code{"numeric"}: The duality gap.}
    \item{\code{rdgap}:}{Object of class \code{"numeric"}: The relative duality gap.}
    \item{\code{certp}:}{Object of class \code{"numeric"}: The
      certificate of primal infeasibility.} 
    \item{\code{certd}:}{Object of class \code{"numeric"}: The
      certificate of dual infeasibility.} 
    \item{\code{pslack}:}{Object of class \code{"numeric"}: The
      smallest primal slack.}    
    \item{\code{dslack}:}{Object of class \code{"numeric"}: The
      smallest dual slack.}    
    \item{\code{status}:}{Object of class \code{"character"}: Whether
      convergence has been achieved \code{"optimal"}, or \code{"unknown"}.}    
    \item{\code{niter}:}{Object of class \code{"integer"}: Count of
      iterations.}
  }
}

\section{Extends}{
Class \code{"\linkS4class{PDV}"} directly.
}

\seealso{
  \code{\linkS4class{PDV}}
}

\keyword{classes}

