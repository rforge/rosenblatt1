\name{importBetaMriImages}
\alias{importBetaMriImages}
\title{Import medical image file}
\usage{
  importBetaMriImages(files, ...)
}
\arguments{
  \item{files}{A list of file names with full path.}

  \item{...}{Arguments passed to
  \code{\link{newMriImageFromFile}}}
}
\value{
  A list containing an \code{\link[tractor.base]{MriImage}} class object
  for each file in the input. 
}
\description{
  Function importing a directory with scans in all formats (except DICOM). 
  The function is essentially a convenience wrapper around \code{\link{newMriImageFromFile}}.
  It receives a list of medical imaging files (in the formats supported by
  \code{\link{newMriImageFromFile}}) and return a list of
  \code{\link{MriImage}} objects-- one for each file
  supplied. 
}
\examples{
\dontrun{
	files<- grep("con.*img", list.files(), value=T)
	files<- sub('\\\\.img', "", files)
	scans<- importBetaMriImages(files, fileType='NIFTI')
	class(scans[[1]])
	}
}
\author{
  Jonathan Rosenblatt
}

