\name{KKTSLV-class}
\Rdversion{1.1}
\docType{class}
\alias{KKTSLV-class}

\title{Class \code{"KKTSLV"}}

\description{
  Class for computing the solution of a KKT-system.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("KKTSLV", ...)} or by
invoking \code{kktSLV1()}.
}

\section{Slots}{
  \describe{
    \item{\code{f}:}{Object of class \code{"function"}: A function for
      solving a KKT-system with arguments \code{x, y, z}.}
    \item{\code{items}:}{Object of class \code{"list"}: Includes
      elements pertinent to the function's body of \code{f}.}
  }
}

\seealso{
  \code{\link{kktSOLVE}}
}

\keyword{classes}

