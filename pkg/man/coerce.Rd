\name{coerce}
\alias{coerce}
\alias{coerce,DEFCP,DEFNL-method}

\title{
  Coercion of convex program objects
}

\description{
  These methods coerce an object of S4-class \code{DEFCP} to an
  S4-object of class \code{DEFNL}. Hereby, the convex nonlinear
  objective is casted in its epigraph form.
}

\usage{
\S4method{coerce}{DEFCP,DEFNL}(from, to = "DEFNL", strict = TRUE)  
}

\arguments{
  \item{from}{An object of S4-class \code{DEFCP}.}  
  \item{to}{A characrter string.}
  \item{strict}{\code{logical}; the returned object must be from the
    target class.}
}

\value{
  An object of S4-class \code{DEFNL}.
}

\seealso{
  \code{\linkS4class{DEFNL}}, \code{\linkS4class{DEFCP}}
}

\keyword{optimize}
