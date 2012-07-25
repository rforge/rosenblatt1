# 
# Author: john.ros@gmail.com
###############################################################################


# Import from analyze format.
beta.arrays<- 
		
# Fit model parameters.
mixture.fit<- mixtureFit(beta.arrays, mixtureControl=mixture.control)

# Plot SPMs  of each parameter.
plot(mixture.fit, ...)

# Compute wilcoxon mask.
significance.mask<- computeMask(mixture.fit, statistic=c("Wilcox", "T"), FDRlevel= )

# Export SPM & mask:

# Overlay prevalence on anatomy.
anatomy<- 
anatomyOverlayPlot(mixtureFit=mixture.fit, anatomy=anatomy , outputFile=)



