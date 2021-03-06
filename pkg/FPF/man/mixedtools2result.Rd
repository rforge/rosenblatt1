\name{mixedtools2result}
\alias{mixedtools2result}
\title{Map mixtools output to simple vector}
\usage{
  mixedtools2result(mixedtools.output, result, model)
}
\arguments{
  \item{mixedtools.output}{Output of
  \code{\link[mixtools]{normalmixEM}}. List of class \emph{mixEM}.}

  \item{result}{The vector to be returned with
  \code{\link[mixtools]{normalmixEM}} output}

  \item{model}{Character vector specifying the type of
  model fitted with \code{\link[mixtools]{normalmixEM}}.}
}
\value{
  Vector of voxel-wise parameter estimates. Utility
  function for \code{\link{pointWiseMixtureFit}}
}
\description{
  Map \code{\link[mixtools]{normalmixEM}} package output to simple vector.
}
\author{
  Jonathan Rosenblatt
}

