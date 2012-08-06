# 
# Author: johnros
###############################################################################

#require(mixtools)
#require(tractor.base)
#
#setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
#(files<- grep("con.*img", list.files(), value=T))
#(files<- sub('\\.img', "", files))
#scans<- importBetaMriImages(files, fileType='NIFTI')
#
#.data<- na.omit(scans[20,20,20,])
#
#temp.fit<- normalmixEM(.data, k=3, mean.constr = c(0,0,NA) )
#names(temp.fit)
#ls.str(temp.fit)
#temp.fit$loglik
#
#temp.fit2<- normalmixEM(.data, k=2, mean.constr = c(0,0) )
#temp.fit2$loglik




##################################################################


mixedtools2result<- function(mixedtools.output, result, model){
	
	# Preparing output for **full** model:
	if(model=="full"){
		component.order<- c( order(mixedtools.output$sigma[c(1,2)]) , 3)
		result[c("p1.1","p2.1","p3.1")] <- mixedtools.output$lambda[component.order] 
		result[c("A.1","B.1","C.1")] <- mixedtools.output$sigma[component.order]^2
		result["mu.1"]<- mixedtools.output$mu[component.order][3]				
	}
	
	
	
	# Preparing output for **full** model:
	if(model=="null"){
		component.order<- order(mixedtools.output$sigma[c(1,2)]) 
		result[c("p1.0","p2.0")] <- mixedtools.output$lambda[component.order] 
		result[c("A.0","B.0")] <- mixedtools.output$sigma[component.order]^2	
	}
		
	
	
	# Preparing output for **full** model:
	return(result)
}
















# Workhorse for fitting mixture models:
pointWiseMixtureFit<- function(beta.vector, fit.control, progress, ...){
	# Verify intput:
	stopifnot(is.numeric(beta.vector))
	
	
	# Initialize output:
	if(!missing(progress)) { setTxtProgressBar(progress, getTxtProgressBar(progress)+1)	}
	
	output.names<-c(
			"initial.p1", "initial.p2", "initial.p3", "initial.mu", "initial.A", "initial.B", "initial.C",
			"p1.1", "p2.1", "p3.1",
			"mu.1",
			"A.1", "B.1", "C.1", 
			"likelihood.1",
			"p1.0","p2.0",
			"A.0","B.0",
			"likelihood.0",
			"n"
	)
	output<- rep(NA, length.out= length(output.names))
	names(output)<- output.names
	n<- sum(!is.na(beta.vector))
	output[['n']]<- n
	if(n < fit.control$minObservations) return(output)
	clean.beta.vector<- na.omit(beta.vector)
	
	
	
	## TODO: B) Add initialization tricks.
	
	# Fit 3 component mixture
	try({
				capture.output(temp.result.three.components <- pointWise3MixtureFit(clean.beta.vector, fit.control))
				output<- mixedtools2result(temp.result.three.components, output, model="full")
				
				
				
				
				
				p3.bound<- p3Bound(p1=output[['p1.1']], mu=output[['mu.1']], A=output[['A.1']], B=output[['B.1']], n=output[['n']], fit.control=fit.control)
				if(output[['p3.1']] < p3.bound) {
					# Note: present implementation forces constraint. Next- optimize with it.
					# TODO: B) Imlpement constrained optimizawtion
					output[['p3.1']]<- p3.bound
				} 
				
				
				
				
				# Fit 2 component mixture
				capture.output(temp.result.two.components<- pointWise2MixtureFit(clean.beta.vector, fit.control, temp.result.three.components))
				output<- mixedtools2result(temp.result.two.components, output, model="null")
				
			}) # End of "try" statement
	
	
	
	
	# Finiliazing:
	return(output)
}
## Testing:
#.data<- scans[20,20,20,]
#(output<- pointWiseMixtureFit(.data, fit.control = generateMixtureControl()))












# Function for pointwise fitting **full** model:
pointWise3MixtureFit<- function(beta.vector, fit.control){
	
	temp.result<- normalmixEM(beta.vector, k=3, mean.constr = c(0,0,NA), verb=FALSE )
	
	return(temp.result)
		
}
## Testing:
#.data<- na.omit(scans[20,20,20,])
#output.3<- pointWise3MixtureFit(.data, fit.control= generateMixtureControl())










# Function for pointwise fitting null model:
pointWise2MixtureFit<- function(beta.vector, fit.control, three.component.fit){

	temp.result<- normalmixEM(beta.vector, k=2, mean.constr = c(0,0), verb=FALSE, lambda=three.component.fit$lambda , mu=three.component.fit$mu, sigma=three.component.fit$sigma)
	
	return(temp.result)
		
}
## Testing:
#.data<- na.omit(scans[20,20,20,])
#output.3<- pointWise3MixtureFit(.data, fit.control= generateMixtureControl())
#(output.2<- pointWise2MixtureFit(.data, fit.control= generateMixtureControl(), three.component.fit = output.3))






