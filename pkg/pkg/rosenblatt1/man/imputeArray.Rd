\name{imputeArray}
\alias{imputeArray}
\title{Utility function}
\usage{
  imputeArray(beta.array, min.neighbours = 20)
}
\arguments{
  \item{beta.array}{A four dimensional array of coefficients. }

  \item{min.neighbours}{Minimum number of valid neighboring observations for the imputation to occur}
}
\value{
  An array with imputed values.

  The function will look for location with more than 20
  neighbours and impute the missing value with their 0.1
  trimmed mean.
}
\description{
  Impute an array with first degree neighbors. Typically called by \code{\link{brainMixtureFit}}.
}
\author{
  Jonathan Rosenblatt
}

