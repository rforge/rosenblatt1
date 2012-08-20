\name{pointWiseMixtureFit}
\alias{pointWiseMixtureFit}
\alias{pointWise2MixtureFit}
\alias{pointWise3MixtureFit}
\title{Pointwise fitting of group mixture}
\usage{
  pointWiseMixtureFit(beta.vector, fit.control, progress)  
  
    pointWise3MixtureFit(beta.vector, fit.control)
    
        pointWise2MixtureFit(beta.vector, fit.control,   three.component.fit) 
}
\arguments{
  \item{beta.vector}{Numeric vector of estimated contrasts.
  One per subject.}

  \item{fit.control}{List of control parameters generated
  by \code{\link{generateMixtureControl}}}

  \item{progress}{A connection object to a progress bar.}
  
  \item{three.component.fit}{An initial fit done by \code{\link{pointWise3MixtureFit}}}
}
\value{
  Numeric vector with estimates of mixture parameters.
}
\description{
  \code{pointWiseMixtureFit} returns the point wise fitting of group mixture.\cr
  \code{pointWise3MixtureFit} returns a \code{mixEM} class object of the \emph{three} component mixture fit.\cr
  \code{pointWise3MixtureFit} returns a \code{mixEM} class object of the \emph{two} component mixture fit.\cr  
}
\author{
  Jonathan Rosenblatt
}

