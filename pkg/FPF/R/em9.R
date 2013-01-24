# TODO: B) Test solver on more data. 
# 
# Author: johnros
###############################################################################

#rm(list=ls())
#require(rosenblatt1)
#dyn.load('/home/johnros/workspace/MixtureRandomEffects/trunk/likelihood.so')
#require(tractor.base)



initializeOutput<- function(){
	output.names<- c(
			#"initial.p1", "initial.p2", "initial.p3", "initial.mu", "initial.A", "initial.B", "initial.C",
			"p1.1", "p2.1", "p3.1",
			"mu.1",
			"A.1", "B.1", "C.1", 
			"likelihood.1",
			"p3.bound",
			"p1.0","p2.0",
			"A.0","B.0",
			"likelihood.0",
			"n"
	)
	output<- rep(NA, length.out= length(output.names))
	names(output)<- output.names
	return(output)
}






rmixednorm<-function(p1, p2, p3, mu, A, B, C, draws){
	group.count<- rmultinom(1, draws, prob=c(p1,p2,p3))
	mus<- c(0,0,mu)
	sds<- sqrt(c(A,B,C))
	result<- c(rnorm(group.count[1], mus[1], sds[1]),rnorm(group.count[2], mus[2], sds[2]) ,rnorm(group.count[3], mus[3], sds[3]) )	
	return(result)
} 
## Testing:
#rmixednorm(0.1, 0.2, 0.7, 1, 1, 2, 1, 10)












generateRandomParams<- function(random.starts, moments){
	arbitrary<- generateArbitraryParams(moments)
	
	result<- with(as.list(arbitrary),{			
				mu.random <-  rnorm(random.starts, mean=0, sd=1)
				A.random <-  abs(A * rexp(n=random.starts, rate=1))
				B.random <-  abs(B * rexp(n=random.starts, rate=1)) 
				C.random <-  abs(C * rexp(n=random.starts, rate=1)) 
				p3.random <-  runif(random.starts) 
				p1.random <-  runif(n=random.starts, min=0 , max=1-p3.random) 
				p2.random <-  1 - p3.random - p1.random
				cbind(p1=p1.random, p2=p2.random, p3=p3.random, mu=mu.random, A=A.random, B=B.random, C=C.random)
			})	
	return(result)
}
## Testing:
#generateRandomParams(10, c(m1=1,m2=1,m3=1,m4=0))






generateArbitraryParams<- function(moments){
	result<-with(as.list(moments),{ 
				p1<- 0.2 
				p2<- 0.2 
				p3<- 1-p1-p2 
				mu<- m1/p3 
				A<- m2/p1 
				B<- m2/p2 
				C<- m2/p3
				return(c(p1=p1,p2=p2,p3=p3,mu=mu,A=A,B=B,C=C))
			})
	return(result)
}
## Testing:
#generateArbitraryParams(c(m1=1, m2=2, m3=1, m4=0))


checkParams<- function(params, n, fit.control){
	ok<- FALSE
	if(isTRUE(all(sapply(params, is.numeric)))){
		ok<- isTRUE(with(as.list(params),{					
					p1 >= 0 && p1<=1 && p2 >= 0 && p2<=1 && p3 >= 0 && p3<=1 &&	A > 0 && B>0 && C>0 && 
							p3 <= p3Bound(mu, A, B, C, n, fit.control) 
				}))		
	}	
	return(ok)
}



generateHybridParams<- function(resolution, moments, fit.control, n){
	m1<- moments[['m1']]
	m2<- moments[['m2']]
	m3<- moments[['m3']]
	m4<- moments[['m4']]
	
	warn.state<- getOption('warn')
	options(warn=-1)
	
	p1<- rep(seq(0, 1, length=resolution), each=resolution)				
	p3<- runif(n =  length(p1), max=1-p1)
	p2<- 1-p3-p1
	
	mu<- m1/p3
	
	C<- (m3-m1^3/p3^2)/(3*m1)		
	A<- (-m1* p1* p3^2* (2* m1^3 - 3* m1* m2* p3 + m3* p3^2) + sqrt(m1^2* p1* p2* p3^3* (2 *m1^6* (p1 + p2 - 2* p3) +
									12* m1^4* m2* p3^2 + 3* m1^2* (-3* m2^2 + m4* (p1 + p2))* p3^3 +
									6* m1* m2* m3* p3^4 - 4* m1^3* m3* p3^2* (p1 + p2 + p3) -
									m3^2* p3^4* (p1 + p2 + p3))))/(3* m1^2* p1* (p1 + p2)* p3^3)
	
	B<- -(m1* p2* p3^2* (2* m1^3 - 3* m1* m2* p3 + m3* p3^2) + sqrt(m1^2* p1* p2* p3^3* (2* m1^6* (p1 + p2 - 2* p3) + 
									12* m1^4* m2* p3^2 + 3* m1^2* (-3* m2^2 + m4* (p1 + p2))* p3^3 + 
									6* m1* m2* m3* p3^4 - 4* m1^3* m3* p3^2* (p1 + p2 + p3) - 
									m3^2* p3^4* (p1 + p2 + p3))))/(3* m1^2* p2 *(p1 + p2)* p3^3)	
	
	C<- pmax(A, C, na.rm=TRUE)
	
	
	temp.result<- cbind(p1=p1, p2=p2, p3=p3, A=A, B=B, C=C, mu=mu)
	check.rows<- apply(temp.result, 1, checkParams, n, fit.control)
	result<- temp.result[check.rows,]
	
	
	options(warn=warn.state)
	return(result)
}
## Testing:
#generateHybridParams(resolution = 10, c(m1=1, m2=2, m3=1, m4=0) , generateMixtureControl(), n=100)
## FIXME: Not working?





# replaces likelihoodem9.c

likelihood.c<- function(p1, p2, p3, mu, A, B, C, beta.vector,...){
	.result<- NA
	
	# Integrity checks:
	if(!(sum(c(p1,p2,p3),na.rm=TRUE)>0)) return(.result)
	
	replace.p1<- 0
	replace.A<- 99 
	replace.p2<- 0
	replace.B<- 99
	replace.p3<- 0
	replace.C<- 99
	replace.mu<- 0
	
	# Components: 1
	if(!is.na(p1)&& !is.na(A)&& is.numeric(p1) && is.numeric(A) && p1<=1 && p1>=0 && A>0){
		replace.p1<- p1
		replace.A<- A
	}
	# Components: 2
	if(!is.na(p2)&& !is.na(B)&& is.numeric(p2) && is.numeric(B) && p2<=1 && p2>=0 && B>0){
		replace.p2<- p2
		replace.B<- B
	}
	# Components: 3
	if(!is.na(p3)&& !is.na(C) &&!is.na(mu)&& is.numeric(p3) && is.numeric(C) && is.numeric(mu) && p3<=1 && p3>0 && C>0){
		replace.p3<- p3
		replace.C<- C
		replace.mu<- mu
	}
	
	.result<- .C("likelihood",
			p1=as.numeric(replace.p1),
			p2=as.numeric(replace.p2),
			p3=as.numeric(replace.p3),
			mu=as.numeric(replace.mu),
			A=as.numeric(replace.A),
			B=as.numeric(replace.B),
			C=as.numeric(replace.C),
			rawData=as.numeric(beta.vector),
			n=length(beta.vector),
			result=as.numeric(0) )
	return(.result$result)
}
## Testing:
#do.call(likelihood.c, list(p1=0.1, p2=0.2, p3=0.7, mu=1, A=1, B=1, C=1, beta.vector=rnorm(100)))
#test.params<- c(p1=0.1, p2=0.2, p3=0.7, mu=1, A=1, B=1, C=1)
#do.call(likelihood.c, c(test.params, list(beta.vector=rnorm(100))) )







makeStartValues<- function(arbitrary.params, random.params, hybrid.params, random.starts, beta.vector, fit.control){
	# Combine random and hybrid starting points:
	n<- length(beta.vector)
	temp.start.values<- rbind(arbitrary.params, random.params, hybrid.params)
	ok.rows<- apply(temp.start.values, 1, checkParams, n, fit.control)
	start.values<- temp.start.values[ok.rows,]
	
	# Now keep only the best results:
	if(nrow(start.values) > random.starts){
		start.likelihoods <- as.numeric(apply(start.values, 1, function(x) {					
							try(do.call(likelihood.c, c(x, list(beta.vector=beta.vector))))									
						}))		
		
		best.ind<- rank(start.likelihoods) > (length(start.likelihoods) - random.starts)
		start.values<- start.values[best.ind, ]		
	}	
	return(start.values)
}
## Testing:
#makeStartValues(
#		arbitrary.params = generateArbitraryParams(c(m1=1,m2=1,m3=1,m4=0)),
#		random.params = generateRandomParams(10, c(m1=1,m2=1,m3=1,m4=0)),
#		hybrid.params = generateHybridParams(resolution = 3, c(m1=1, m2=2, m3=1, m4=0), n=10 ),
#		5,
#		beta.vector,
#		generateMixtureControl())











arrangeVariances<- function(params){
	result<- within(as.list(params),{				
				if (A > B){
					temp<- p1
					p1<- p2
					p2<- temp
					temp<- B
					B<- A
					A<- temp
				}				
			})
	return(unlist(result))	
}






















#replaces initializeEM9:
initialize3MixtureFitFast<- function(beta.vector, fit.control){
	comparePrecision<- fit.control$roundTolerance
	thresh<- fit.control$numericThresh
	iteration.limit<- fit.control$iteration.limit
	random.starts<- fit.control$RandomStarts
	resolution<- fit.control$Resolution
	n<- length(beta.vector)
	
	#Calculating sample moments
	m1<- mean(beta.vector)
	m2<- mean(beta.vector^2)
	m3<- mean(beta.vector^3)	
	m4<- mean(beta.vector^4)
	empirical.moments<- c(m1=m1, m2=m2, m3=m3, m4=m4)
	
	# Setting initial arbitrary values:
	arbitrary.params<- generateArbitraryParams(empirical.moments)
	
	
	# But also try other starting points...
	random.params<-generateRandomParams(random.starts=random.starts, empirical.moments)	
	
	
	
	# Try grid with moment estimators:
	hybrid.params<- generateHybridParams(resolution=resolution, moments=empirical.moments, fit.control, n=n)
	
	
	# Combine all starting candidates and remove "illegal" options:
	start.values<- makeStartValues(arbitrary.params, random.params, hybrid.params, random.starts = random.starts, beta.vector = beta.vector, fit.control)	
	
	# Choose best starting value	
	temp.result<- selectBestInitializtion(start.values=start.values, beta.vector, fit.control=fit.control)	
	
	# Assure B is the larger variance component:
	result<- arrangeVariances(temp.result)
		
	return(result)
}
## Testing:
#beta.vector<- rmixednorm(0.1, 0.2, 0.7, 1, 1, 2, 1, 10)
#initialize3MixtureFitFast(beta.vector, generateMixtureControl())











# *Non* Null Responsability function:
responsibility3<- function(params, beta.vector, fit.control){
	
	
	outer.result<- with(as.list(params), {
				expo1<- -0.5 * beta.vector^2 / A
				expo1[is.na(expo1)]<- fit.control$numericThresh
				phi1<- 1/sqrt(A) * exp(expo1)
				
				expo2<- -0.5 * beta.vector^2 / B
				expo2[is.na(expo2)]<- fit.control$numericThresh
				phi2<- 1/sqrt(B) * exp(expo2)
				
				expo3<- -0.5 * (beta.vector-mu)^2 / C
				expo3[is.na(expo3)]<- fit.control$numericThresh
				phi3<- 1/sqrt(C) * exp(expo3)
				
				# Only third component:
				if(isTRUE(round(p3, fit.control$roundTolerance)==1)){
					total<- p3*phi3
					inner.result<- list(resp1=0/total, resp2=0/total, resp3=total/total)
					return (inner.result)			
				}
				# Only first component:
				if(isTRUE(round(p1, fit.control$roundTolerance)==1)){
					total<- p1*phi1
					inner.result<- list(resp1=total/total, resp2=0/total, resp3=0/total)
					return (inner.result)			
				}
				# Only second component:
				if(isTRUE(round(p2, fit.control$roundTolerance)==1)){
					total<- p2*phi2
					inner.result<- list(resp1=0/total, resp2=total/total, resp3=0/total)
					return(inner.result)			
				}
				# Only first two components:
				if(isTRUE(round(p3, fit.control$roundTolerance)==0)){
					total<- p1*phi1 + p2*phi2 
					inner.result<- list(resp1=p1*phi1/total, resp2=p2*phi2/total, resp3=0/total)	
					return (inner.result)			
				}
				
				
				total<- p1*phi1 + p2*phi2 + p3*phi3
				inner.result<- list(resp1=p1*phi1/total, resp2=p2*phi2/total, resp3=p3*phi3/total)
				
				return (inner.result)	
			}
	
	)
	return(outer.result)
}
## Testing:
#beta.vector<- rmixednorm(0.1, 0.2, 0.7, 1, 1, 2, 1, 100)
#responsibility3(c(p1=0.1, p2=0.2, p3=0.7, mu=1, A=1, B=1, C=1), beta.vector, generateMixtureControl())
#responsibility3(c(p1=0.8, p2=0.2, p3=0, mu=1, A=1, B=1, C=1), beta.vector, generateMixtureControl())
#table(unlist(responsibility3(c(p1=0.8, p2=0.2, p3=0, mu=1, A=1, B=1, C=1), beta.vector, generateMixtureControl()))<0)










logit<- function(p) log(p/(1-p))

inv.logit<- function(x) {
	if(x > log(.Machine$double.xmax)) return(1)
	exp(x) / (1+exp(x))
}
## Testing:
#inv.logit(logit(0.9999999999999999995))
#inv.logit(500)
#names(options())












# Maximization step:
m.step3<-function(resps, constrained, beta.vector, fit.control){
	## Initialization:
	variance.bound<- fit.control$variance.bound
	mStep.iteration.limit<- fit.control$mStep.iteration.limit
	result<- c(p1=NA, p2=NA, p3=NA, mu=NA, A=NA, B=NA, C=NA)
	n<- length(beta.vector)
	
	resps.total<- sum(unlist(resps))
	T0.1<- sum(resps[['resp1']])
	T0.2<- sum(resps[['resp2']])
	T0.3<- sum(resps[['resp3']])
	T1.1<- resps[['resp1']]%*%beta.vector
	T1.2<- resps[['resp2']]%*%beta.vector
	T1.3<- resps[['resp3']]%*%beta.vector
	T2.1<- resps[['resp1']]%*%beta.vector^2
	T2.2<- resps[['resp2']]%*%beta.vector^2
	T2.3<- resps[['resp3']]%*%beta.vector^2	
	
	
	# Unconstrained estimates:	
	result[['p1']]<- T0.1 / resps.total
	result[['p2']]<- T0.2 / resps.total
	result[['p3']]<- 1-result[['p1']]-result[['p2']]			
	result[['mu']]<- T1.3 / T0.3  
	result[['A']]<- min(T2.1 / T0.1, variance.bound)
	result[['B']]<- min(T2.2 / T0.2, variance.bound)
	result[['C']]<- min(max(result[['A']], T2.3 / T0.3 - result[['mu']]^2), variance.bound)			
	
	
	
	# Deal with constrained optimization (just improve, no need to find optimum)
	if(constrained){
		
		
		constrained.target.function<- function(args){			
			p1<- inv.logit(args[['logit.p1']]) 
			mu<- args[['mu']] 
			A<- exp(args[['log.A']]) 
			B<- exp(args[['log.B']]) 
			C<- exp(args[['log.C']])
						
			# Here is the difference with past constraint functions:
			expo<- exponentialConstraint(mu,A,B,C,n)
			
			result<- -1e15
			
			# Make sure constraint on p1 is enforced during initialization and iteration
			if(is.nan(args[['logit.p1']]) || is.infinite(args[['logit.p1']]) || log(p1) > expo    ) return(result)			
			
			value<- c(-0.5*n*log(2*pi)-
							0.5*log(A)*T0.1 + 	log(p1)*T0.1 -	0.5/A * T2.1 -
							0.5*log(B)*T0.2 + 	log(exp(expo)-p1)*T0.2 - 	0.5/B * T2.2 -
							0.5*log(C)*T0.3 + 	log(1-exp(expo))*T0.3 - 0.5/C*(T2.3 + mu^2*T0.3 - 2*mu*T1.3))			
			
			if( isTRUE(!is.na(value) && !is.infinite(value)) ) result<- value	
			
			return(result)
		}
										
		# Initialization values for optimization:
		init.par<- with(as.list(result), list(logit.p1=logit(p1), mu=mu, log.A=log(A), log.B=log(B), log.C=log(C))   )
		
		## Actual optimization:		
		## TODO: A) Deal initialization values that return -Inf!
		optim.result<- optim(
				par=init.par,		
				fn=constrained.target.function,				
				control=list(fnscale=-1, maxit=mStep.iteration.limit ), method='Nelder-Mead')		
		
		
		result[['mu']]<- optim.result$par[['mu']]
		result[['A']]<- exp(optim.result$par[['log.A']])
		result[['B']]<- exp(optim.result$par[['log.B']])
		result[['C']]<- exp(optim.result$par[['log.C']])
		result[['p1']]<- inv.logit(optim.result$par[['logit.p1']])
		result[['p3']]<- p3Bound(mu=result[['mu']], A=result[['A']], B=result[['B']], C=result[['C']], n=n, fit.control = fit.control)
		result[['p2']]<- 1 - result[['p3']] - result[['p2']]				
	}
	
	return(result)
}
## Testing:
#beta.vector<- rmixednorm(0.5,0.3,0.2,1,1,2,1,10000)
#m.step3(
#		responsibility3(c(p1=0.5, p2=0.3, p3=0.2, mu=1, A=0.2, B=2, C=2.0001), beta.vector, generateMixtureControl()),
#		constrained=TRUE, 
#		beta.vector, 
#		generateMixtureControl())
	



















# replaces iterate.em9
# returns a vector of parameter values OR NAs if it did not converge
iterate3MixtureFitFast<- function(initial.params, beta.vector, iteration.limit, constrained, fit.control){
	# Assumes beta.vactor has no NAs
	
	## Initialization:
	n<- length(beta.vector)
	change<- 1
	iterations<- 0L
	percision.tolerance<- fit.control$iterationTolerance
	roundTolerance<- fit.control$roundTolerance
	mStep.iteration.limit <- fit.control$mStep.iteration.limit
	new.params<- arrangeVariances(initial.params)		
	
	while( change > percision.tolerance && iterations < iteration.limit){
		
		resps<- responsibility3(new.params, beta.vector, fit.control)
		
		if(all(round(resps[['resp3']], roundTolerance)==0)){
			# There are only two components:			
			return(iterate2MixtureFitFast(new.params, beta.vector, iteration.limit, fit.control))			
			# Note: this solution does not allow to estimate 3 componets once resps3 have vanished.
		}
		
		## Calculate new parameter values:
		old.params<- new.params		
		
		new.params<- m.step3(resps, constrained, beta.vector, fit.control)
		
		likeOld<- do.call(likelihood.c, c(old.params, list(beta.vector=beta.vector)) )					
		likeNew<- do.call(likelihood.c, c(new.params, list(beta.vector=beta.vector)) )
		change<- abs(likeNew-likeOld)
		
		iterations<- iterations+1L			
	}		
	
	return(c(new.params, likelihood=likeNew))	
}
## Testing:
#beta.vector<- rmixednorm(0.7,0.3,0,0,1,2,1,100)
#iterate3MixtureFitFast(c(p1=0.7, p2=0.2, p3=0.1, mu=1, A=1, B=2, C=1), beta.vector, 20, constrained = TRUE, generateMixtureControl())
#load('debugging.Rdata')
#iterate3MixtureFitFast(initial, beta.vector, 20, constrained = TRUE, generateMixtureControl())






check.location<- function(params, beta.vector, fit.control){
	result<- rep(NA, length(params))
	names(result)<- names(params)	
	try(
			result<- iterate3MixtureFitFast(
					initial.params=params, 
					beta.vector = beta.vector, 
					fit.control = fit.control, 
					iteration.limit = 5, # number of iterations before choosing initialization point 
					constrained = FALSE)
			, silent=TRUE)
	return(result)	
}
## Testing:
check.location(c(0.1, 0.3, 0.6,  1, 1, 2, 1), beta.vector, fit.control)







check.likelihood<- function(params, beta.vector){
	result<- NA
	if(isTRUE(all(sapply(params, is.numeric)))){
		result<- do.call(likelihood.c, c(params, list(beta.vector=beta.vector)))
	}			
	return(result)
}
## Testing:










selectBestInitializtion<- function(start.values, beta.vector, fit.control){
	# From each starting value
	# iterate a little
	# keep end values
	
	end.values<- t(apply(start.values, 1, check.location, beta.vector, fit.control))
	if(ncol(end.values)==1L) { end.values<- t(end.values) }
	
	# compute likelihood
	# return most likelly values
	
	temp.check.likelihood<- function(params) check.likelihood(params, beta.vector)
	likelihoods <- apply(end.values, 1, temp.check.likelihood)	
	likelihoods<- as.numeric(likelihoods)
	row.of.max<- which.max(likelihoods)
	
	return(end.values[row.of.max,])
	
}
## Testing:
#beta.vector<- FPF:::rmixednorm(0.1,0.2,0.7,1,1,2,1,100)
selectBestInitializtion(t(as.matrix(c(p1=0.1,p2=0.3,p3=0.7, mu=1,A=1,B=2,C=1), byrow=TRUE)), beta.vector, generateMixtureControl() )
















# Replace pointWise3MixtureFit()
pointWise3MixtureFitFast<- function(beta.vector, fit.control){
	# beta.vector is assumed to have no NAs.
	
	n<- length(beta.vector)		
	
	
	## Initializing:
	initial<- initialize3MixtureFitFast(beta.vector, fit.control = fit.control)		


	## First-- try to solve unconstrained:
	temp.result<-iterate3MixtureFitFast(
			initial.params=initial, 
			beta.vector=beta.vector,			 
			iteration.limit=fit.control$iteration.limit, 
			constrained = FALSE,
			fit.control = fit.control)
	
	
	
	# If the bound on p3 is effective:
	p3.bound<- p3Bound(mu=temp.result[['mu']], A=temp.result[['A']], B=temp.result[['B']], C=temp.result[['C']],  n=n, fit.control=fit.control)
	if(temp.result[['p3']] > p3.bound ){		
		
		temp.result<- iterate3MixtureFitFast(
				initial.params=initial,		
				beta.vector=beta.vector,				 
				iteration.limit=fit.control$iteration.limit, 
				constrained = TRUE,
				fit.control = fit.control)
	}	
	
	
	return(temp.result)		
}
## Testing:
## FIXME: Constrained estimation not working.
#hist(beta.vector<- rmixednorm(0.5, 0.2, 0.3, 2, 0.5, 0.2, 0.2, 1000))
#round(pointWise3MixtureFitFast(beta.vector, generateMixtureControl()),2)






# like responsability3:
responsibility2<- function(params, beta.vector, fit.control){
	# assumed params if given in the format of a 3(!) component mixture
	stopifnot(params[['p3']]==0)
	result<- responsibility3(params, beta.vector, fit.control)
	return(result)
}
## Testing:
#responsibility2(c(p1=0.5,p2=0.5,p3=0,mu=1,A=1,B=2,C=1), beta.vector, generateMixtureControl())
#.resps<- responsibility2(initialize2with3(pointWise3MixtureFitFast(beta.vector, generateMixtureControl())), beta.vector, generateMixtureControl())



	
	
	
	

# like m.step3:	
m.step2<- function(resps, beta.vector, fit.control){
	result<- m.step3(resps, constrained = FALSE, beta.vector, fit.control)
	if(round(result[['p3']],fit.control$roundTolerance)==0)  result[['p3']]<- 0
	else stop("p3 under the null is not 0 (returned by m.step3)")
	return(result)
}
## Testing:
#beta.vector<- rmixednorm(0.1,0.2,0.7,1,1,2,1,50)
#.resps<- responsibility2(initialize2with3(pointWise3MixtureFitFast(beta.vector, generateMixtureControl())), beta.vector, generateMixtureControl())






initialize2with3<- function(mixture3result){
	stopifnot(is.numeric(mixture3result))
	# returns a vector of parameters with p3=0
	outer.result<- with(as.list(mixture3result),{
				p3<- 0
				new.p1<- p1 / (p1+p2)
				new.p2<- 1-new.p1
				inner.result<- c(p1=new.p1, p2=new.p2, p3=p3, mu=NaN, A=A, B=B, C=NaN)
				return(inner.result)
			})
	result <- arrangeVariances(outer.result)
	return(unlist(result))
	 
}
## Testing:
#beta.vector<- rmixednorm(0.6,0.1,0.3,1,1,2,1,100)
#initialize2with3(pointWise3MixtureFitFast(beta.vector, generateMixtureControl()))







# replaces iterate.em9
iterate2MixtureFitFast<- function(mixture3result, beta.vector, iteration.limit, fit.control){
	# Assumes beta.vactor has no NAs
	# Is initialized with the output of pointWise3MixtureFitFast()
	
	## Initialization:
	n<- length(beta.vector)
	change<- 1
	iterations<- 0
	percision.tolerance<- fit.control$iterationTolerance
	roundTolerance<- fit.control$roundTolerance
	mStep.iteration.limit <- fit.control$mStep.iteration.limit
	new.params<- initialize2with3(mixture3result)		
	
	while( change > percision.tolerance && iterations < iteration.limit){
		
		resps<- responsibility2(new.params, beta.vector, fit.control)
		
		## Calculate new parameter values:
		old.params<- new.params		
		
		new.params<- m.step2(resps, beta.vector, fit.control)
		
		likeOld<- do.call(likelihood.c, c(old.params, list(beta.vector=beta.vector)) )					
		likeNew<- do.call(likelihood.c, c(new.params, list(beta.vector=beta.vector)) )
		change<- abs(likeNew-likeOld)
		
		iterations<- iterations+1			
	}		
	
	return(c(new.params, likelihood=likeNew))	
}
## Testing:
#beta.vector<- rmixednorm(0.1,0.2,0.7,1,1,2,1,100)
#iterate2MixtureFitFast(c(p1=0.7, p2=0.2, p3=0.1, mu=1, A=1, B=2, C=1), beta.vector, 20, generateMixtureControl())















# Replace pointWise2MixtureFit()
pointWise2MixtureFitFast<- function(beta.vector, fit.control, mixture3result){
	# beta.vector is assumed to have no NAs.	
	n<- length(beta.vector)		
	
	
	temp.result<-iterate2MixtureFitFast(
			mixture3result=mixture3result, 
			beta.vector=beta.vector,			 
			iteration.limit=fit.control$null.iteration.limit,
			fit.control = fit.control
			)	
	
	return(temp.result)	
}
## Testing:
#beta.vector<- rmixednorm(0.6,0.4,0,1,1,2,1,5000)
#pointWise2MixtureFitFast(beta.vector, generateMixtureControl(), pointWise3MixtureFitFast(beta.vector, generateMixtureControl()))











# Reaplce mixedtools2result()
fit2result<- function(mixture.fit, result, model=c("full","null")){
	
	# Preparing output for **full** model:
	if(model=="full"){
		result[c("p1.1","p2.1","p3.1","mu.1","A.1","B.1","C.1","likelihood.1")] <- mixture.fit[c("p1","p2","p3","mu","A","B","C","likelihood")]			
	}
	
	
	
	# Preparing output for **null** model:
	if(model=="null"){
		result[c("p1.0","p2.0","A.0","B.0","likelihood.0")] <- mixture.fit[c("p1","p2","A","B","likelihood")]	
	}
	
	return(result)
}
## Testing:











# Fit voxel-wise null and full mixture 
pointWiseMixtureFitFast<- function(beta.vector, fit.control, progress=NULL){
	# Verify intput:
	stopifnot(is.numeric(beta.vector))
	
	
	# Initialize output:
	if(!missing(progress)) { setTxtProgressBar(progress, getTxtProgressBar(progress)+1)	}
	
	output<-initializeOutput()
	
	n<- sum(!is.na(beta.vector))
	output[['n']]<- n
	if(n < fit.control$minObservations) return(output)
	clean.beta.vector<- na.omit(beta.vector)
	
	
	
	try({
				# Fit 3 component mixture (constraint implemented within pointWise3MixtureFitFast)
				temp.result.three.components <- pointWise3MixtureFitFast(clean.beta.vector, fit.control)
				output<- fit2result(temp.result.three.components, output, model="full")				
				
				
				# Fit 2 component mixture
				temp.result.two.components<- pointWise2MixtureFitFast(clean.beta.vector, fit.control, temp.result.three.components)
				output<- fit2result(temp.result.two.components, output, model="null")
				
			}, silent=TRUE) 	
	
	# Finiliazing:
	return(output)
}
## Testing:
#beta.vector<- rosenblatt1:::rmixednorm(0.4,0.4,0.2,1,1,1,1, 50)
#beta.vector<- rmixednorm(0.4,0.4,0.2,1,1,1,1, 50)
#pointWiseMixtureFitFast(beta.vector, generateMixtureControl())
#data(VinkData)
#pointWiseMixtureFitFast(MriImage2Array(scans)[20,20,20,], generateMixtureControl())







	
	
brainMixtureFitArray<- function(beta.array, fit.control){
	
	dims<-dim(beta.array)
	cat('This may take several minutes. Why not load the fortunes package and enjoy those quotes?\n')
	progress.bar<- txtProgressBar(min=0, max=prod(dims[-4]), style=1)
	warn <- options(warn = 2)
	
	pointwise<- function(beta.vector) {
		unlist(pointWiseMixtureFitFast(beta.vector = beta.vector, fit.control = fit.control, progress=progress.bar))
	}
	first.fit <- apply(beta.array, c(1,2,3), pointwise)
	close(progress.bar)
	
	
	
	# TODO: B) Add imputation and smoothing before re-estimation + check which solution has highest likelihood
	## Smooth estimates:
	# impute array before smoothing
	# smooth estimates
	
	
	smoothed.fit<- first.fit
	
	## Second fit:
	second.fit <- smoothed.fit	
	
	
	## Impute non-convergence
	imputed.fit<-  wrapImputeArray(second.fit)
	
	return(imputed.fit)
}
## Testing:
#load('/home/johnros/workspace/MixtureRandomEffects/pkg/rosenblatt1/data/VinkData.RData')
#fit.array<- brainMixtureFitArray(MriImage2Array(scans)[20:22,20:22,20:22,], generateMixtureControl())	




	
















brainMixtureFit<- function(MRImage.list, fit.control= generateMixtureControl()){
	
	## Verify input:
	stopifnot(is.list(MRImage.list) && class(MRImage.list[[1]])=="MriImage" )	
	
	
	## Logging:
	#log.file.name<- tempfile("Mixture_Log_")
	#log.con<- file(log.file.name, open="wt")	
	#cat("Log written to", log.file.name,"\n")
	warn<- options('warn'=0)
	
	
	## Initialize:
	beta.array<- MriImage2Array(MRImage.list)	
			
	## Fit	
	imputed.fit<- brainMixtureFitArray(beta.array, fit.control)
	
	## Finilizing
	dim.names<- dimnames(imputed.fit)
	fitted.list<- list()
	meta<- newMriImageMetadataFromTemplate(MRImage.list[[1]]$getMetadata(), datatype=getDataTypeByNiftiCode(16))
	for(param in dim.names[[1]]){
		fitted.list[[param]]<- newMriImageWithData(imputed.fit[param,,,], meta)
	}
	
	
	## Exit:
	options(warn)
	message('Done')
	return(fitted.list)	
}
## Testing:
#require(rosenblatt1)

#test.brain.fit<- brainMixtureFit(scans, generateMixtureControl())
#	
#x11()
#createSliceGraphic(test.brain.fit[["p3.1"]], z=26)
#image(test.brain.fit[["p3.1"]]$getData()[,,26])
##save(test.brain.fit, file="/home/johnros/workspace/MixtureRandomEffects/pkg/rosenblatt1/data/VinkDataFit.RData")
##load(file="/home/johnros/workspace/Mixture Random Effects/tmp/VinkData.RData")
#lapply(test.brain.fit, function(x) x$getData()[20,20,20])


