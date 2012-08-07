# 
# Author: john.ros@gmail.com
###############################################################################
require(tractor.base)
require(rosenblatt1)

# Import SPM{betas} for all subjects:
setwd('~/Projects/MRI/Data/Vink/data_vink_2010/cons_no_smooth')
(files<- grep("con.*img", list.files(), value=T))
(files<- sub('\\.img', "", files))
beta.arrays<- importBetaMriImages(files, fileType='NIFTI')
image(beta.arrays[,,20,21])



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



