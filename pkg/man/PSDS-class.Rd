\name{PSDS-class}
\Rdversion{1.1}
\docType{class}
\alias{PSDS-class}

\title{Class \code{"PSDS"}}

\description{
  Nesterov-Todd scalings and scaled Lagrange-multipliers.
}

\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PSDS",
  ...)} or by using \code{nnos()}.
}

\section{Slots}{
  \describe{
    \item{\code{W}:}{Object of class \code{"list"} with elements
      \code{r} and \code{rti} for the computation of the Nesterov-Todd
      scaling and its inverse and the scaled variable \code{lambda}.}
  }
}

\examples{
showClass("PSDS")
}
