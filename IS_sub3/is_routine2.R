
#ismoments(target,sigma.cav,as.numeric(mu.cav),numnit,tol,abctable,rej=F)
ismoments = function(target,sigma,mu,is.iter,isss) 
#mu is the mean vector from the cavity distribution - we treat it as the prior mean for log(theta) and log(rho)
#sigma is the covariance matrix for the cavity distribution - we treat it as the prior covarance for log(theta) and log(rho)

#NB this code is inefficient because we can actually re-use an initially simulated ABC table, without resimulating
#for each locus, but this just gets things going.
{
	nparams = length(mu) # = 3 in this case with log(theta), log(R), log(tf), 
	
	trans.sigma = sigma
	diag(trans.sigma) = 2*diag(trans.sigma)
	params = rmnorm(is.iter, mean = mu, varcov = trans.sigma) #simulate from widened cavity distribution (acts as a prior)
	                                                      #use a widened distribution to improve stability
	                                                      #correct it below in simwt
	simwt = dmnorm(params,mean=mu,varcov=sigma)/dmnorm(params,mean=mu,varcov=trans.sigma)
	write(t(params),"paramfile",ncol=nparams)
	


	use_ms_tf_scaling=T #see runnit.r for explanation
	cm1 = paste("./is_moments",target,as.integer(is.iter),as.integer(isss),"paramfile",as.integer(use_ms_tf_scaling))
	wvec = system(cm1,intern=T)
	wvec = as.numeric(wvec)
	wvec = simwt*wvec #apply correction
	ess = sum(wvec)^2/sum(wvec^2) #the effective sample size 
	pss = ess/is.iter

	if(length(params[,1])!=length(wvec)){
		print(paste("length params ",length(params[,1])))
		print(paste("length wvec ",length(wvec)))
		stop("error in length of wvec")
	}
	
	if(ess <= nparams + 2){
		print(paste("ess problem for",target,": returning prior"))
		est.mu = mu
		est.cov = sigma
		est.precmat = chol2inv(chol(sigma))
		est.sigmainv.mu = as.numeric(est.mu%*%est.precmat)
 		outlist = list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,
 						sigma = est.cov,ess = ess,pss=pss,params=params,wvec=wvec)	
 		return(outlist)	
		
	}
	
	est.mu = sweep(params,1,wvec,"*")
	est.mu = apply(est.mu,2,sum)/sum(wvec)
	
	
	
	p.cent = sweep(params,2,est.mu,"-")	#!!!I forgot this earlier - took hours to find!!!
	p.cent.left = sweep(p.cent,1,wvec,"*")
	est.cov = (t(p.cent.left) %*% p.cent)/sum(wvec)
	#calc above should be equivalent to est.cov = (t(p.cent) %*% diag(wvec) %*% p.cent)/sum(wvec)
	#need to do it this way for efficiency...
    est.cov = est.cov*ess/(ess - 1) 

	est.precmat = chol2inv(chol(est.cov)) #find the inverse
	est.precmat = (ess-nparams - 2)/(ess-1)*est.precmat #this gives an unbiased estimate of precision matrix - again using wss in 
	                                                    #place of actual sample size...
	
 	est.sigmainv.mu = as.numeric(est.mu%*%est.precmat) 
 	print(target)
	print("est.mu")
	print(est.mu)
	print("est.cov")
	print(est.cov)
	print("est.precmat")
	print(est.precmat)

 	list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,sigma = est.cov,ess = ess,pss=pss,params=params,wvec=wvec)		


}
