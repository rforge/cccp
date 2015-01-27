\name{gets}
\alias{gets}
\alias{gets,PDV-method}
\alias{gets,CPS-method}

\title{
  Extractor for slot \code{s}
}

\description{
  Returns the value of the list-variable(s) \code{s}.
}

\usage{
\S4method{gets}{PDV}(object)  
\S4method{gets}{CPS}(object)  
}

\arguments{
  \item{object}{An object of either S4-class \code{PDV} or \code{CPS}.}
}

\value{
  \code{list}, the value of the list-variable(s) \code{s}.
}

\keyword{optimize}
