
# Author: Jonathan Rosenblatt john.ros@gmail.com
###############################################################################
##' Fit mixture models to fMRI group studies.
##' 
##' \tabular{ll}{
##' Package: \tab rosenblatt1\cr
##' Type: \tab Package\cr
##' Version: \tab 0.5\cr
##' Date: \tab 2012-08-19\cr
##' License: \tab GPL (>= 2)\cr
##' }
##' 
##' Mixture model fitting for random effect in fMRI group studies.
##'
##' @name rosenblatt1-package
##' @aliases rosenblatt1
##' @docType package
##' @title Fit mixture models to fMRI group studies
##' @include tractor.base
##' @include mixtools
##' @author Jonathan Rosenblatt 
#NA






#
##' Function imoprting a directory with scans in all formats except DICOM:
##' 
##' @param files A list of file names with full path. 
##' @param ... Arguments passed to \code{\link{newMriImageFromFile}}
##' @return A list containin an \code{\link{MriImage}} class object for each file in the input.
##' 
##' The function is a convinience wrapper around \code{\link{newMriImageFromFile}}. It recieves a list of medical imaging files (in the formats supported by \code{\link{newMriImageFromFile}}) and return a list of \code{\link{MriImage}} objects-- one for each file supplied.
##' 
##' @author Jonathan Rosenblatt 
##' @examples 
##' files<- grep("con.*img", list.files(), value=T)
##' files<- sub('\\.img', "", files)
##' scans<- importBetaMriImages(files, fileType='NIFTI')
##' class(scans[[1]])
##' @export
importBetaMriImages<- function(files, ...){
	# Verify input:
	#~~~
	
	
	# Import:
	require(tractor.base)
	subject.arrays<- list()
	.test.dimentions<- newMriImageFromFile(fileName =files[1], ...)$getDimensions()
	for(i in seq_along(files)){
		subject.arrays[[i]]<- newMriImageFromFile(files[i], ...)
		# Test all scans have the same dimensions:
		if(any(.test.dimentions!= subject.arrays[[i]]$getDimensions())) { stop('Arrays with different resolutions')}
	}
	
	return(subject.arrays)
	
	
}
## examples:
#setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
#(files<- grep("con.*img", list.files(), value=T))
#(files<- sub('\\.img', "", files))
#scans<- importBetaMriImages(files, fileType='NIFTI')
#class(scans[[1]])=="MriImage"
#save(scans, file= '/home/johnros/workspace/Mixture Random Effects/truck/rosenblatt1/data/VinkData.Rdata')























#
##' Generates a list of control arguments.
##' @param numericThresh Numeric value.  
##' @param minObservations Minimal number of observation to attempt fitting.
##' @return List with control parameters for \code{\link{brainMixtureFit}}
##' 
##' @author Jonathan Rosenblatt 
##' @export
generateMixtureControl<- function(
		numericThresh= -18, 
		minObservations= 15, # minimum number of observations for estimation
		iteration.limit=100,
		null.iteration.limit=1000,
		roundTolerance= 15,
		RandomStarts= 5,
		Resolution= 20,
		mStep.iteration.limit=5,
		variance.bound=10,
		iterationTolerance=1e-4,
    seed=333
){	
	return(list(
					numericThresh= numericThresh, 	
					minObservations= minObservations,
					iteration.limit=iteration.limit,
					null.iteration.limit=null.iteration.limit,
					roundTolerance=roundTolerance,
					RandomStarts=RandomStarts,
					Resolution=Resolution,
					mStep.iteration.limit=mStep.iteration.limit,
					variance.bound=variance.bound,
					iterationTolerance=iterationTolerance,
          seed=seed
			))	
}

	
	

























#
##' Generate matrix of neighbouring indexes for brain image imputation
##' 
##' @param i 
##' @param j 
##' @param k 
##' @return Matrix of neighbouring voxel indexes.
##' 
##' @author Jonathan Rosenblatt 
neighbourhoodMatrix<- function(i,j,k){
	as.matrix(expand.grid(
					c(i-1,i,i+1), 
					c(j-1,j,j+1), 
					c(k-1,k,k+1))[-14,])
}
## Testing:
#rosenblatt1:::neighbourhoodMatrix(10,20,30)










# Deprecated by imputeArray2 
imputeArray <- function(beta.array, min.neighbours=20) {
	
	result<- beta.array
	dims<- dim(result)
	
	for(i in 2:(dims[1]-1)){
		for(j in 2:(dims[2]-1)){
			for(k in 2:(dims[3]-1)){
				neighbour.matrix<- neighbourhoodMatrix(i,j,k)
				if( is.na(beta.array[i,j,k]) && sum(!is.na(beta.array[neighbour.matrix])) > min.neighbours ){
					result[i,j,k]<- mean(beta.array[neighbour.matrix], trim=0.1, na.rm=TRUE)					
				}
			}
		}
	}
	
	return(result)
}
## Testing:
#test.array<- array(rnorm(10*10*10), dim=c(10,10,10))
#rosenblatt1:::imputeArray(test.array)






# Impute an array with first degree neighbours (Shlomi's C implementation):
#dyn.load("/home/johnros/workspace/MixtureRandomEffects/trunk/average_NA.so")

imputeArray2 <- function(beta.array, min.neighbours=20) {	
	dims<- dim(beta.array)	
	
		result <- .Call("fill_NA_voxels_with_neig_mean",
				x=as.double(beta.array),
				dim.x=as.integer(dims[1]),
				dim.y=as.integer(dims[2]),
				dim.z=as.integer(dims[3]))
		
		return (array(data=result,dim=dims))
}
## Testing:
#test.array<- array(rnorm(10*10*10), dim=c(10,10,10))
#imputeArray2(test.array)




















#
##' Impute a selected parameters in a brain array
##' @param brain.mixture.fit.object 
##' @param param.to.impute 
##' @return An array of parameter estimates imputed using first degree neighbourhood.
##' 
##' @author Jonathan Rosenblatt 
wrapImputeArray<- function(brain.mixture.fit.object){
	## Initializing:
	dims<- dim(brain.mixture.fit.object)
	dim.names<- dimnames(brain.mixture.fit.object)		
	result<- array(NA, dim=dims, dimnames = dim.names)
	
	
	for (param in dim.names[[1]]){
		result[param,,,]<- imputeArray2(brain.mixture.fit.object[param,,,])
			}
			
	return(result)
	
}
## Testing:
#my.array<- array(rnorm(10^4), dim=c(x=10,y=10,z=10,j=10), dimnames = list(LETTERS[1:10], NULL, NULL, NULL ))
#dimnames(my.array)
#wrapImputeArray(my.array)










p3Bound<-function(p1, p2, mu, A, B, C, n, fit.control){
	stopifnot(is.numeric(c(p1, p2, mu, A, B, C, n)))
	
	result<- 0
	expo<- fit.control$numericThresh
	
	try({
				expo<- -( mu^2 / (2 * (p1*A+p2*B)) )
			}, silent=TRUE)
	
	# Detection version: 
	 if(!is.na(expo)) result<-   exp(expo) / sqrt(n)
	
	# Clasification version ("estimation" in Donoho & Jin):
	#if(!is.na(expo)) result<-   exp(expo) 
	
	return(result)
}
## Test:
#p3Bound(0.2, 0.6, 1,1,2,1.01,67, generateMixtureControl())







# Returns TRUE if p3 in **permissible** zone.
checkBound <- function(p1, p2, p3,mu, A, B, C, n, fit.control){
	p3 >= p3Bound(p1, p2, mu, A, B, C, n, fit.control)
}
## Testing:
#checkBound(0.2, 0.2, 0.6, 1, 1,2,1,60,generateMixtureControl())
#checkBound(p1= 0.9, p2= 0.09, p3= 0.01, mu= 1, A=1, B=2, C=1, n=60,generateMixtureControl())






#
##' Converts a list of MRIimage objects to a single four dimentional array.
##' @param MRImage.list List of Mriimage objects. All assumed to have same resolutino and representing different subjects.
##' @return A numeric four dimentional array with space as first three dimentions and subject as the fourth.
##' 
##' @author Jonathan Rosenblatt 
##' @export
##' @examples
##' data(VinkData)
##' MriImage2Array(scans)

MriImage2Array<- function(MRImage.list){
	
	
	## Initialize: 
	single.scan.dim<- MRImage.list[[1]]$getDimensions()
	number.of.sbjects<- length(MRImage.list)
	dims<-c(single.scan.dim, number.of.sbjects)
	
	
	beta.array<- array(0, dim=dims, dimnames = c('x','y','z', 'subjects'))
	for(i in seq_along(MRImage.list)){
		beta.array[, , , i]<- MRImage.list[[i]]$getData()		
	}
	
	# Finalize:
	return(beta.array)
	
}
## Testing:
#setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
#(files<- grep("con.*img", list.files(), value=T))
#(files<- sub('\\.img', "", files))
#scans<- importBetaMriImages(files, fileType='NIFTI')
#MriImage2Array(scans)



























#
##' Export all estimated parameters to imaging files
##' @param mixture.fit.object Output of \code{\link{brainMixtureFit}}
##' @param format The medical imaging format to export. See details.
##' @param file.heading Heading of output files.
##' @return Nothing. Called for it's side effects.
##' 
##' This is a convenience wrapper around \code{\link{writeMriImageToFile}}.
##'  
##' It exports all the SPM outputed by \code{\link{brainMixtureFit}} and saves them in any of the formats supported by \code{\link{writeMriImageToFile}. 
##' @author Jonathan Rosenblatt 
##' @export
##' @examples
##' #exportFitAsMedicalImage(mixture.fit.object = test.brain.fit, file.heading = "export", format = "test")
exportFitAsMedicalImage<- function(mixture.fit.object, format, file.heading){
	## Verify input:
	stopifnot(is.character(file.heading))
	stopifnot(is.character(format))
	stopifnot(all(sapply(mixture.fit.object, class)=="MriImage"))
	
	
	## Initialize:
	params<- names(mixture.fit.object)	
	for(param in params){
		file.name<- paste(file.heading, "_", param, sep="")
		writeMriImageToFile(mixture.fit.object[[param]], fileName = file.name, fileType=format, overwrite = FALSE)
	}
		
	## Finalize:
	cat("Files save in ", getwd())
}
## Testing:
#setwd("/home/johnros/workspace/Mixture Random Effects/tmp/")
#exportFitAsMedicalImage(mixture.fit.object = test.brain.fit, file.heading = "export", format = "dfsdfs")















#
##' Computes a rejecton group mask given beta maps.
##' @param MRImage.list List of Mriimage objects. All assumed to have same resolutino and representing different subjects.
##' @param test.statistic Which group test statistic to use. At present "T" or "Wilcoxon" (default) are supported.
##' @param FDR.level Level of FDR multiplicity control. 
##' @param fit.control List of control parameters generated by \code{\link{generateMixtureControl}} 
##' @return An \link{MriImage} class object corresponding to the binary mask of rejected locations.
##' 
##' @author Jonathan Rosenblatt 
##' @export
##' @examples
##' data(VinkDataFit) 
##' test.brain.mask<- computeMask(MRImage.list=scans, method = "T", fit.control = generateMixtureControl())
##' createSliceGraphic(test.brain.mask, z=30)

computeMask<- function(MRImage.list, test.statistic="Wilcoxon", FDR.level=0.1, fit.control=generateMixtureControl()){
	## Verify input:
	stopifnot(is.list(MRImage.list) && class(MRImage.list[[1]])=="MriImage" ) 
	stopifnot(is.numeric(FDR.level) && FDR.level<1 && FDR.level>0)
	
	
	
	## Initialize:
	beta.array<- MriImage2Array(MRImage.list)
	dims<- dim(beta.array)	
	min.n<- fit.control$minObservations
	
	
	
	# Create map of p-values according to the method
	if(test.statistic=="T"){
		pvals.array<- apply(beta.array, c(1,2,3), function(x){
					result<- NA
					if(sum(is.finite(x)) < min.n) return(result)
					try(result <- t.test(x)[['p.value']] )
					return(result)
				}	) 
	}
	
	
	
	else if(test.statistic=="Wilcoxon"){
		pvals.array<- apply(beta.array, c(1,2,3), function(x){
					result<- NA
					if(sum(is.finite(x)) < min.n) return(result)
					try(result<- wilcox.test(x, mu=0, alternative="two.sided", exact=FALSE)$p.value)
					return(result)
				}	)
	}
	
		
	else{
		message("The method supplied is not valid.")
	}
	
	
	
	## Compute FDR adjusted p-values
	output<- array(p.adjust(pvals.array, method = 'BH')<= FDR.level, dim=dims[1:3])
	
	
	## Finialize:
	mri.mask<- newMriImageWithData(output, MRImage.list[[1]])
	return(mri.mask)
	
	
}
## Testing:
# require(tractor.base)
# setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
# files<- grep("con.*img", list.files(), value=T)
# files<- sub('\\.img', "", files)
# scans<- importBetaMriImages(files, fileType='NIFTI')
# test.brain.mask<- computeMask(MRImage.list=scans, fit.control = generateMixtureControl())
# createSliceGraphic(test.brain.mask, z=30)


























