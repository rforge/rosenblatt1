# 
# Author: Jonathan Rosenblatt john.ros@gmail.com
###############################################################################



########## Testing tractor.base ########
#test<-do.call(newMriImageFromFile, list(
#				fileName ='/home/johnros/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth/ruis290510_1_con_0001', 
#				fileType = 'ANALYZE'))
#test$getData()
#test$getOrigin()
#
#test2<-do.call(newMriImageFromFile, list(
#				fileName ='/home/johnros/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth/ruis290510_1_con_0001' 
#				))
#test2$getData()
#test2$getOrigin()
################################################################



# Function imoprting a directory with scans in all formats except DICOM:
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
		if(any(.test.dimentions!= subject.arrays[[i]]$getDimensions())) { error('Arrays with different resolutions')}
	}
	
	
	#Arranging data in one array:
	(all.scans.dimensions<-c(.test.dimentions, length(files))) 
	all.scans<- array(0, dim=all.scans.dimensions, dimnames = c('x','y','z', 'subjects'))
	for(i in seq_along(files)){
		all.scans[, , , i]<- subject.arrays[[i]]$getData()		
	}
		
	return(all.scans)
	
	
}
## Example:
#setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
#(files<- grep("con.*img", list.files(), value=T))
#(files<- sub('\\.img', "", files))
#scans<- importBetaMriImages(files, fileType='NIFTI')
#image(scans[,,20,21])





## TODO: B) Create import function for DICOM format.
# Function imoprting a directory with scans in DICOM format:








# Creates an empty vector with givens alot names:
initizalizeOutput<- function(the.names){
	stopifnot(is.character(the.names))	
	output<- rep(NA, times=length(the.names))
	names(output)<- the.names
	return(output)
}






generateMixtureControl<- function(
		numericThresh= -18, # ???
		iteration.limit= 100, # number of EM iterations for fitting the full model
		finalizeIterations= 1000,  # number of iteratinos to fit the null model
		roundTolerance= 15,	# rounding precision for numerican comparisons (faster than all.equal() )
		minObservations= 15, # minimum number of observations for estimation
		iterationTolerance= 1e-4, # minimal change in parameter estimates (or likelihood?) that terminates iteration. 
		Resolution= 10, # resolutino of grid of initialization values for p1,p2
		NullResolution= 7, # resolutino of initialization grid for null estimation
		RandomStarts= 5,  # number of random initializatino values to try
		initiateIterations=5,  # number of iteration on all initializatino points before picking the best and finish iterating
		NullRandomStarts= 10, # number of iteration on all *null* initializatino points before picking the best and finish iterating
		step.shrink=0.75, # the shrinkage of the p3 estimate from the constrained value when the constraint holds.
		recursionDepth= 0, # ?
		recursionLimit= 1, # ?
		varianceBound= 30 # bound on the value of component
){	
	return(list(
					numericThresh= numericThresh, iteration.limit= iteration.limit,	
					finalizeIterations= finalizeIterations,roundTolerance= roundTolerance,	
					minObservations= minObservations, iterationTolerance= iterationTolerance, 
					Resolution= Resolution, NullResolution= NullResolution, 
					RandomStarts= RandomStarts, initiateIterations= initiateIterations, 
					NullRandomStarts= NullRandomStarts, step.shrink = step.shrink, 
					recursionDepth= recursionDepth, recursionLimit= recursionLimit
			))	
}

	
	









## Main fitting function (if mixedtools is not used)
#pointWiseMixtureFit<- function(beta.vector, output, null.fit, progress.bar, fit.control, ...){
#	# Initialize
#	
#	# Remove missing observations
#	
#	
#	# If left with enough observations-- fit mixture.
#	# If p3 constaint violated: refit while forcing constraint. 
#	# Note: unclear how to fit while constraining weights...
#	
#}




# Utility function to create the index matrix of neighobours:
neighbourhoodMatrix<- function(i,j,k){
	as.matrix(expand.grid(c(i-1,i,i+1), c(j-1,j,j+1), c(k-1,k,k+1))[-14,])
}




#' Impute an array with first degree neighbours
#' @param beta.array 
#' @param min.neighbours 
#' @returnType 
#' @return 
#' 
#' @author johnros
#' @export
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












#' Impute a selected parameters in a brain array
#' @param brain.mixture.fit.object 
#' @param param.to.impute 
#' @returnType 
#' @return 
#' 
#' @author johnros
#' @export
wrapImputeArray<- function(brain.mixture.fit.object){
	dims<- dim(brain.mixture.fit.object)
	dim.names<- dimnames(brain.mixture.fit.object)
		
	result<- array(NA, dim=dims, dimnames = dim.names)
	
	for (param in dim.names[[4]]){
		result[,,,param]<- imputeArray(brain.mixture.fit.object[,,,param])
			}
			
	return(result)
	
}
## Testing:
#my.array<- array(rnorm(10^4), dim=c(x=10,y=10,z=10,j=10), dimnames = list(NULL, NULL, NULL, LETTERS[1:10]))
#dimnames(my.array)
#wrapImputeArray(my.array)













#' Checks if constraint on p3 is met?
#' @param p1 
#' @param mu 
#' @param A 
#' @param B 
#' @param n 
#' @param fit.control 
#' @returnType 
#' @return 
#' 
#' @author johnros
#' @export
p3Bound<-function(p1, mu, A, B, n, fit.control){
	result<- 1
	stopifnot(is.numeric(c(p1,A, B, mu)))
	expo<- fit.control$numericThresh
	try(expo<- -abs(mu) * sqrt(n / (A + B)), silent=TRUE)
	if(!is.na(expo)) result<- 1 - exp(expo)			
	return(result)
}
## Test:
#p3.bound(0.2,0.5,3,1,2,67,EM9.control)




















#' Fits a full and null mixture models to a SPM{beta}
#' @param beta.array  Four diemntional array of fitted betas
#' @param fit.control List of control parameters generated by \link{\code{generateMixtureControl}} 
#' @param ... 
#' @returnType 
#' @return 
#' 
#' @author johnros
#' @export
brainMixtureFit<- function(beta.array, fit.control= generateMixtureControl(), ...){
	## TODO: B) Finish mixture fitting function.
	## TODO: A) Use S4-MRIimage or S3 class object to allow export 
	
	# Verify input:
	stopifnot(length(dim(beta.array))==4L)
	
	
		
	# Initialize:
	dims<- dim(beta.array)
	message('This may take several minutes. Why not load the fortunes package and enjoy those quotes?')
	progress.bar<- txtProgressBar(min=0, max=prod(dims[-4]), style=3)
	warn <- options(warn = 2)
	
	
	
	
	
	## Fit- first step
	pointwise<- function(beta.vector) {
		unlist(pointWiseMixtureFit(beta.vector = beta.vector, fit.control = fit.control, progress=progress.bar))
	}
	
	
	first.fit <- apply(beta.array, c(1,2,3), pointwise) 
		
	
	
	
	## TODO: B) Add imputation and smoothing before re-estimation + check which solution has highest likelihood
	
	## Smooth estimates:
	# impute array before smoothing
	# smooth estimates	
	smoothed.fit<- first.fit	
	
	
	
	## Second fit:
	second.fit<- smoothed.fit
	
	
	
	
	## Impute non-convergence
	imputed.fit<-  wrapImputeArray(second.fit)
	
	


	## Finilizing
	
	close(progress.bar)
	options(warn)
	message('Done')
	
	class(imputed.fit)<- "mixtureSPM"	
	return(imputed.fit)
	
	
}
## Testing:
#require(tractor.base)
#rm(list=ls())
#setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
#files<- grep("con.*img", list.files(), value=T)
#files<- sub('\\.img', "", files)
#beta.arrays<- importBetaMriImages(files, fileType='NIFTI')
#test.brain.fit<- brainMixtureFit(beta.arrays[25:30,28:30,28:30,], fit.control = generateMixtureControl())
#dim(test.brain.fit)
#dimnames(test.brain.fit)
#table(is.na(test.brain.fit))
#image(test.brain.fit["p3.1",,,30])
#save(test.brain.fit, file="/home/johnros/workspace/Mixture Random Effects/tmp/test.brain.fit.Rdata")
#load(file="/home/johnros/workspace/Mixture Random Effects/tmp/test.brain.fit.Rdata")


















# Function to export estimated parameters
exportArrayAsMedicalImage<- function(mixture.fit.object, format, file){
	## Verify input:
	stopifnot(is.character(file))
	
	
	## Initialize:
	require(tractor.base)
	
	
	param.list<- c("p1.1")
	for (param in )
	
	
	
	newMriImageFromTemplate(mixture.fit.object)
	
	
	
	# save MRImage object
	
	
	
	
	## Finalize:
}















# Function to compute rejection mask using either group-t or Wilcoxon test statistics
computeMask<- function(beta.array, method=c("T","Wilcoxon"), FDR.level=0.1, fit.control){
	# Verify input:
	stopifnot(is.numeric(beta.array) && isTRUE(length(dim(beta.array))==4L))
	stopifnot(is.numeric(FDR.level) && FDR.level<1 && FDR.level>0)
	
	
	
	# Initialize:
	dims<- dim(beta.array)	
	min.n<- fit.control$minObservations
	
	
	
	# Create map of p-values according to the method
	if(method=="T"){
		pvals.array<- apply(beta.array, c(1,2,3), function(x){
					result<- NA
					if(sum(is.finite(x)) < min.n) return(result)
					try(result <- t.test(x)[['p.value']] )
					return(result)
				}	) 
	}
	
	
	
	else if(method=="Wilcoxon"){
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
	
	
	
	# Compute FDR adjusted p-values
	output<- array(p.adjust(pvals.array, method = 'BH')<= FDR.level, dim=dims[1:3])
	
	
	# Finialize:
	return(output)	
	
	
}
## Testing:
#require(tractor.base)
#setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
#files<- grep("con.*img", list.files(), value=T)
#files<- sub('\\.img', "", files)
#beta.arrays<- importBetaMriImages(files, fileType='NIFTI')
#test.brain.mask<- computeMask(beta.array=beta.arrays, method = "T", fit.control = generateMixtureControl())
#image(test.brain.mask[,,30])
















# Function to plot selected parameters at a given slice:
SPMSlicePlot<- function(mixture.fit.array, parameter, slice){
	
}







# Function to plot selected parameters overlayn on anatomy given a mask: 
AnatomyOverayPlot<- function(mixture.fit.array, anatomy.array, parameter, slice){
	
}



