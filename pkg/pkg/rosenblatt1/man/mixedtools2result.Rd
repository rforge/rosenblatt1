\name{mixedtools2result}
\alias{mixedtools2result}
\title{Map mixtools output to simple vector}
\usage{
  mixedtools2result(mixedtools.output, result, model)
}
\arguments{
  \item{mixedtools.output}{Output of
  \code{\link{normalmixEM}}. List of class \emph{mixEM}.}

  \item{result}{The vector to be returned with
  \code{\link{normalmixEM}} output}

  \item{model}{Character vector specifying the type of
  model fitted with \code{\link{normalmixEM}}.}
}
\value{
  Vector of voxel-wise parameter estimates. Utility
  function for \code{\link{pointWiseMixtureFit}}
}
\description{
  Map mixtools output to simple vector
}
\author{
  Jonathan Rosenblatt
}

