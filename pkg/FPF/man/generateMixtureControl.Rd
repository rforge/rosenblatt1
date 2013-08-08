\name{generateMixtureControl}
\alias{generateMixtureControl}

\title{
Generates control object for fitting.
}
\description{
A constructor function for the list of arguments needed by \code{\link{brainMixtureFit}}
}
\usage{
generateMixtureControl(numericThresh = -18, minObservations = 15, 
  iteration.limit = 100, null.iteration.limit = 1000, roundTolerance = 15, 
  RandomStarts = 5, Resolution = 20, mStep.iteration.limit = 5, 
  variance.bound = 10, iterationTolerance = 1e-04, seed=333
)
}

\arguments{
  \item{numericThresh}{The percision used for numerical comparisons}
  
  \item{minObservations}{Minimal number of observations.}
  
  \item{iteration.limit}{Number of iterations to estimate full model}
  
  \item{null.iteration.limit}{Number of iteration to estimate null model}
  
  \item{roundTolerance}{Rounding tolerance for numerical comparisons}
  
  \item{RandomStarts}{Number of random starting values}
  
  \item{Resolution}{Number of starting grid values}
  
  \item{mStep.iteration.limit}{Number of iteration in the \emph{constrained} m-step}
  
  \item{variance.bound}{Upper bound on the values of variance parameters}
  
  \item{iterationTolerance}{Magnitude of change required before convergence}
  
  \item{seed}{Fixed the seed of the number sequence for initialization}
}
\value{Returns a list with the above parameters.}

\author{
Jonathan Rosenblatt <john.ros.work@gmail.com>
}

