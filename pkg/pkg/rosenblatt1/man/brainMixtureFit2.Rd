\name{brainMixtureFit}
\alias{brainMixtureFit}
\title{Fits a mixture models to a group SPM{beta}}
\usage{
  brainMixtureFit(MRImage.list,
    fit.control = generateMixtureControl())
}
\arguments{
  \item{MRImage.list}{List of Mriimage objects. All assumed
  to have same resolution and representing different
  subjects.}

  \item{fit.control}{List of control parameters generated
  by \code{\link{generateMixtureControl}}}
}
\value{
  List of MriImage objects. Each holding the SPM of a
  different parameters: 
  \item{p1.1}{Estimate of \emph{full} mixture parameters.}
  \item{p2.1}{Estimate of \emph{full} mixture parameters.} 
  \item{p3.1}{Estimate of \emph{full} mixture parameters.} 
  \item{mu.1}{Estimate of \emph{full} mixture parameters.} 
  \item{A.1}{Estimate of \emph{full} mixture parameters.} 
  \item{B.1}{Estimate of \emph{full} mixture parameters.} 
  \item{C.1}{Estimate of \emph{full} mixture parameters.} 
  \item{likelihood.1}{\emph{Full} maximal likelihood.}
  \item{p3.bound}{Deprecated.} 
  \item{p1.0}{Estimate of \emph{null} mixture parameters.} 
  \item{p2.0}{Estimate of \emph{null} mixture parameters.} 
  \item{A.0}{Estimate of \emph{null} mixture parameters.} 
  \item{B.0}{Estimate of \emph{null} mixture parameters.}
  \item{likelihood.0}{\emph{Null} maximal likelihood}
  \item{n}{Sample size.} 
}
\description{
  Given a list of MriImage objects-- all registered to a
  common template-- the function will fit a mixture model
  at each voxel and return an array of the voxel-wise
  parameter estimates. The mixture being fitted is
  parametrized as: \deqn{p_1 \phi_{0,\sigma_1}+p_2
  \phi_{0,\sigma_2}+p_3 \phi_{\mu,\sigma_3} } where
  \eqn{\phi_{\mu,\sigma}} is the Gaussian PDF.
}
\examples{
data(VinkData)
\dontrun{
	No need to run. 
	Output already in VinkDataFit in order to save time.
	test.brain.fit<- brainMixtureFit(scans, fit.control = generateMixtureControl())
}
data(VinkDataFit)
createSliceGraphic(test.brain.fit[["p3.1"]], z=26)
image(test.brain.fit[["p3.1"]]$getData()[,,26])
lapply(test.brain.fit, function(x) x$getData()[20,20,20])
}
\author{
  Jonathan Rosenblatt
}
\references{
  Revisiting Multi-Subject Random Effects in fMRI
  Rosenblatt J.D., Vink M., Benjamini Y. (In Review)
}

