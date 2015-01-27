\name{getz}
\alias{getz}
\alias{getz,PDV-method}
\alias{getz,CPS-method}

\title{
  Extractor for slot \code{z}
}

\description{
  Returns the value of the list-variable(s) \code{z}.
}

\usage{
\S4method{getz}{PDV}(object)  
\S4method{getz}{CPS}(object)  
}

\arguments{
  \item{object}{An object of either S4-class \code{PDV} or \code{CPS}.}
}

\value{
  \code{list}, the value of the list-variable(s) \code{z}.
}

\keyword{optimize}
