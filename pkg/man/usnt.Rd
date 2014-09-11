\name{usnt}
\alias{usnt}
\alias{usnt,Matrix,NLFS-method}
\alias{usnt,Matrix,NNOS-method}
\alias{usnt,Matrix,SOCS-method}
\alias{usnt,Matrix,PSDS-method}

\title{
  Nesterov-Todd scaling of cone constraints and variables.
}

\description{
  This function applies the Nesterov-Todd scaling to cone constraints or
  cone variables. 
}

\usage{
\S4method{usnt}{Matrix,NLFS}(u, W, trans = FALSE, inv = FALSE)  
\S4method{usnt}{Matrix,NNOS}(u, W, trans = FALSE, inv = FALSE)  
\S4method{usnt}{Matrix,SOCS}(u, W, trans = FALSE, inv = FALSE)  
\S4method{usnt}{Matrix,PSDS}(u, W, trans = FALSE, inv = FALSE)  
}

\arguments{
\item{u}{An object of pertaining either to NNO or SOC or PSD constraints.}
\item{W}{An object of either S4-class \code{NLFS}, or \code{NNOS}, or
  \code{SOCS}, or \code{PSDS}.}
\item{trans}{\code{logical}, whether the transpose of the scaling matrix
  should be used.}
\item{inv}{\code{logical}, whether the inverse of the scaling matrix
  should be used.}
}

\value{
  An S4-object of the corresponding class for cone variables,
  \emph{e.g.}, \code{NLFV}, \code{NNOV}, or \code{SOCV} or
  \code{PSDV}.  
}

\keyword{optimize}
