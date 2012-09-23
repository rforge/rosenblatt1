\name{importBetaMriImages}
\alias{importBetaMriImages}
\title{Function imoprting a directory with scans in all formats except DICOM:}
\usage{
  importBetaMriImages(files, ...)
}
\arguments{
  \item{files}{A list of file names with full path.}

  \item{...}{Arguments passed to
  \code{\link{newMriImageFromFile}}}
}
\value{
  A list containin an \code{\link{MriImage}} class object
  for each file in the input.

  The function is a convinience wrapper around
  \code{\link{newMriImageFromFile}}. It recieves a list of
  medical imaging files (in the formats supported by
  \code{\link{newMriImageFromFile}}) and return a list of
  \code{\link{MriImage}} objects-- one for each file
  supplied.
}
\description{
  Function imoprting a directory with scans in all formats
  except DICOM:
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

