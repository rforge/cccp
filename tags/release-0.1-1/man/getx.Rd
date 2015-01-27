\name{getx}
\alias{getx}
\alias{getx,PDV-method}
\alias{getx,CPS-method}

\title{
  Extractor for slot \code{x}
}

\description{
  Returns the value of the variable \code{x}.
}

\usage{
\S4method{getx}{PDV}(object)  
\S4method{getx}{CPS}(object)  
}

\arguments{
  \item{object}{An object of either S4-class \code{PDV} or \code{CPS}.}
}

\value{
  \code{numeric}, the value of the variable \code{x}.
}

\keyword{optimize}
