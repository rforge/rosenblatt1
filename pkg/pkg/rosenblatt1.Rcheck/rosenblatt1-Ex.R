pkgname <- "rosenblatt1"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rosenblatt1')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("MriImage2Array")
### * MriImage2Array

flush(stderr()); flush(stdout())

### Name: MriImage2Array
### Title: Converts a list of MRIimage objects to a single four dimentional
###   array.
### Aliases: MriImage2Array

### ** Examples

data(VinkData)
MriImage2Array(scans)



cleanEx()
nameEx("brainMixtureFit")
### * brainMixtureFit

flush(stderr()); flush(stdout())

### Name: brainMixtureFit
### Title: Fits a mixture models to a group SPMbeta
### Aliases: brainMixtureFit

### ** Examples

data(VinkData)
## Not run: 
##D 	### No need to run. Output already in VinkDataFit in order to save time.
##D 	test.brain.fit<- brainMixtureFit(scans, fit.control = generateMixtureControl())
## End(Not run)
data(VinkDataFit)
createSliceGraphic(test.brain.fit[["p3.1"]], z=26)
image(test.brain.fit[["p3.1"]]$getData()[,,26])
lapply(test.brain.fit, function(x) x$getData()[20,20,20])



cleanEx()
nameEx("computeMask")
### * computeMask

flush(stderr()); flush(stdout())

### Name: computeMask
### Title: Computes a rejection group mask given beta maps.
### Aliases: computeMask

### ** Examples

data(VinkData)
test.brain.mask<- computeMask(MRImage.list=scans, test.statistic = "T", fit.control=generateMixtureControl())
createSliceGraphic(test.brain.mask, z=30)



cleanEx()
nameEx("exportFitAsMedicalImage")
### * exportFitAsMedicalImage

flush(stderr()); flush(stdout())

### Name: exportFitAsMedicalImage
### Title: Export all estimated parameters to imaging files
### Aliases: exportFitAsMedicalImage

### ** Examples

## Not run: 
##D 	exportFitAsMedicalImage(mixture.fit.object = test.brain.fit, file.heading = "export", format = "test")
## End(Not run)



cleanEx()
nameEx("importBetaMriImages")
### * importBetaMriImages

flush(stderr()); flush(stdout())

### Name: importBetaMriImages
### Title: Function imoprting a directory with scans in all formats except
###   DICOM:
### Aliases: importBetaMriImages

### ** Examples

## Not run: 
##D 	files<- grep("con.*img", list.files(), value=T)
##D 	files<- sub('\\.img', "", files)
##D 	scans<- importBetaMriImages(files, fileType='NIFTI')
##D 	class(scans[[1]])
##D 	
## End(Not run)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
