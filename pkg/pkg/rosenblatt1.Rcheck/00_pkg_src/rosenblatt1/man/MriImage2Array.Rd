\name{MriImage2Array}
\alias{MriImage2Array}
\title{Converts a list of MRIimage objects to a single four dimentional array.}
\usage{
  MriImage2Array(MRImage.list)
}
\arguments{
  \item{MRImage.list}{List of Mriimage objects. All assumed
  to have same resolutino and representing different
  subjects.}
}
\value{
  A numeric four dimentional array with space as first
  three dimentions and subject as the fourth.
}
\description{
  Converts a list of MRIimage objects to a single four
  dimentional array.
}
\examples{
data(VinkData)
MriImage2Array(scans)
}
\author{
  Jonathan Rosenblatt
}

