\name{wrapImputeArray}
\alias{wrapImputeArray}
\title{Utility function}
\usage{
  wrapImputeArray(brain.mixture.fit.object)
}
\arguments{
  \item{brain.mixture.fit.object}{A four dimensional array of mixture parameter estimates.}

}
\value{
  An array of parameter estimates imputed using first
  degree neighborhood.
}
\description{
  Impute a selected parameters in a brain array. To be called by \code{\link{brainMixtureFit}}.
}
\author{
  Jonathan Rosenblatt
}

