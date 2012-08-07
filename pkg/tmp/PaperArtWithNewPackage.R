# 
# Author: johnros
###############################################################################


##### Comparison with my own EM looks good!!!! #####



# Convert MriImage p3 object to p3 array:

source("/home/johnros/workspace/Zero Inflation/R code/PreparePaperArt.R", echo=FALSE, encoding="UTF-8")
setwd("/home/johnros/workspace/Mixture Random Effects/tmp/")
load(file="/home/johnros/workspace/Mixture Random Effects/tmp/VinkData.Rdata")
p3.array<- test.brain.fit$p3.1$getData()
spm.slice<- 26L
heatmap.p3.thresh<- 0.5
p3.thresh.mask<- ifelse(p3.array>heatmap.p3.thresh ,1, NA)
plot.anatomy.overlay(spm=p3.array, anatomy=anatomy2.array, heatmap.mask=p3.thresh.mask, slice= spm.slice, p3.thresh=0, alpha= 0.8, 
		spm.file.out='test.spm', label.file.out = 'test.labels')




#### Demonstrate mixture locations:
source("/home/johnros/workspace/Zero Inflation/R code/PreparePaperArt.R", echo=FALSE, encoding="UTF-8")
ind.array<- p3.array>0.9 & p3.array<1 & apply(scans, c(1,2,3), length)>40

tail(mixture.locations<- arrayInd(which(ind.array), .dim=dim(p3.array)))
dim(mixture.locations)



getParameterFromLocation<- function(MriFit, location){
	i<-mixture.locations[loc,1]
	j<- mixture.locations[loc,2]
	k<- mixture.locations[loc,3]			
	output<- as.list(lapply(MriFit, function(x) x$getData()[i,j,k]))
	return(output)

}


{
	select.locations<- c(200, 1400, 3100,  310)
#	postscript(file='effect.in.select.locations.eps')
	x11()
	par(mfrow=c(2,2))
	for (loc in select.locations){
		#  loc<- 310
		.data<- scans[mixture.locations[loc,1], mixture.locations[loc,2], mixture.locations[loc,3], ]
		hist(.data, freq=FALSE, border='darkgrey', main="", ylim=c(0,1.5), 
				xlab=paste("Location= (",mixture.locations[loc,1],",",mixture.locations[loc,2],",",mixture.locations[loc,1],")", sep="")
		)
		rug(.data)
		
		.parameters<- getParameterFromLocation(test.brain.fit , mixture.locations[loc,])
				
		
		lines(density(x=.data, na.rm=TRUE, bw=0.25), col='darkgrey', ylim=c(0,1.5), xlab="", xlim=c(-3,3),
				main=paste("Prevalence=",round(.parameters$p3.1,2), sep="")	)
		
			
		
				
		
		curve(add=TRUE, .parameters$p1.1*dnorm(x, mean = 0, sd=sqrt(.parameters$A.1)), lty=3)
		curve(add=TRUE, .parameters$p2.1*dnorm(x, mean = 0, sd=sqrt(.parameters$B.1)), lty=3)
		curve(add=TRUE, .parameters$p3.1*dnorm(x, mean = .parameters$mu.1, sd=sqrt(.parameters$C.1)), lty=3)
		mixture<- function(x){
			.parameters$p1.1*dnorm(x, mean = 0, sd=sqrt(.parameters$A.1)) +
					.parameters$p2.1*dnorm(x, mean = 0, sd=sqrt(.parameters$B.1)) +
					.parameters$p3.1*dnorm(x, mean = .parameters$mu.1, sd=sqrt(.parameters$C.1))										
		}
		curve(add=TRUE, mixture)
		abline(v=.parameters$mu.1, lty=2)
		t.statistic<- t.test(.data)$statistic
		title(sub = paste("Group t-Statistic=",round(t.statistic,3)))
	}
		 
}
dev.off()