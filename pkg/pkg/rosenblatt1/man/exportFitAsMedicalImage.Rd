\name{exportFitAsMedicalImage}
\alias{exportFitAsMedicalImage}
\title{Export all estimated parameters to imaging files}
\usage{
  exportFitAsMedicalImage(mixture.fit.object, format,
    file.heading)
}
\arguments{
  \item{mixture.fit.object}{Output of
  \code{\link{brainMixtureFit}}}

  \item{format}{The medical imaging format to export. See
  details.}

  \item{file.heading}{Heading of output files.}
}
\value{
  Nothing. Called for it's side effects.

  This is a convenience wrapper around
  \code{\link{writeMriImageToFile}}.

  It exports all the SPM outputed by
  \code{\link{brainMixtureFit}} and saves them in any of
  the formats supported by
  \code{\link{writeMriImageToFile}}.
}

\description{
  Export all estimated parameters to imaging files
}
\examples{
\dontrun{
	exportFitAsMedicalImage(mixture.fit.object = test.brain.fit, file.heading = "export", format = "test")
}
}
\author{
  Jonathan Rosenblatt
}

