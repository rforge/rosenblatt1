# Author Shlomi Lifshits (lifshits_il@yahoo.com)
# Date: September 2012 

res <- dyn.load("/home/johnros/workspace/MixtureRandomEffects/trunk/average_NA.so")

impute.arr<- function(x)
  {
    res <- .Call("fill_NA_voxels_with_neig_mean",
                 x=as.double(x),
                 dim.x=as.integer(dim(x)[1]),
                 dim.y=as.integer(dim(x)[2]),
                 dim.z=as.integer(dim(x)[3]))
    return (array(data=res,dim=dim(x)))
  }

arr <- array(runif(64),dim=c(4,4,4))
arr[23] <- NA

## BEFORE IMPUTATION
arr

## AFTER IMPUTATION
system.time(impute.arr(arr))

