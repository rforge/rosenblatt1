\name{computeMask}
\alias{computeMask}
\title{Computes rejection group mask}

\usage{
  computeMask(MRImage.list, test.statistic = "Wilcoxon", FDR.level = 0.1, fit.control=generateMixtureControl())
}

\arguments{
  \item{MRImage.list}{List of \link[tractor.base]{MriImage} objects. All assumed
  to have same resolution and representing different
  subjects.}

  \item{test.statistic}{Which group test statistic to use.
  At present "T" or "Wilcoxon" (default) are supported.}

  \item{FDR.level}{Level of FDR multiplicity control.}

  \item{fit.control}{List of control parameters generated
  by \code{\link{generateMixtureControl}}}
}

\value{
  An \link[tractor.base]{MriImage} class object corresponding to the
  binary mask of rejected locations.
}

\description{
  Computes a rejection group mask given beta maps.
}

\examples{
data(VinkData)
test.brain.mask<- computeMask(MRImage.list=scans, test.statistic = "T", fit.control=generateMixtureControl())
createSliceGraphic(test.brain.mask, z=30)
}

\author{
  Jonathan Rosenblatt
}

