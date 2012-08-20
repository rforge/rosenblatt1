\name{generateMixtureControl}
\alias{generateMixtureControl}
\title{Generates a list of control arguments.}
\usage{
  generateMixtureControl(numericThresh = -18,   minObservations = 15)
}
\arguments{
  \item{numericThresh}{Numeric value.}

  \item{minObservations}{Minimal number of observation to
  attempt fitting.}
}
\value{
  List with control parameters for
  \code{\link{brainMixtureFit}}
}
\description{
  Generates a list of control arguments.
}
\author{
  Jonathan Rosenblatt
}

