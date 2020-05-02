# GAUSSIAN AVERAGED EXPECTATION PROPAGATION
library(MASS)
#install.packages(list.files(pattern="tar.gz"), repos = NULL)
library(parallel)
library(mnormt)
library(mcmc)

setwd("/home/mauro/Documents/University/PhDMiniproject/parallelizing")



ismoments = function(target,sigma,mu,is.iter,isss) #numnit=5000 nisamp=20
{
  nparams = length(mu) 
  trans.sigma = sigma
  diag(trans.sigma) = 2*diag(trans.sigma)
  params = rmnorm(is.iter, mean = mu, varcov = trans.sigma) #simulate from widened cavity distribution (acts as a prior)
  #use a widened distribution to improve stability
  #correct it below in simwt
  simwt = dmnorm(params,mean=mu,varcov=sigma)/dmnorm(params,mean=mu,varcov=trans.sigma)
  paramfile <- paste("paramfiles/paramfile_new", ilocus, sep="")
  write(t(params), paramfile,ncol=nparams)
  
  use_ms_tf_scaling=T #see runnit.r for explanation
  cm1 = paste("./is_moments_new",target,as.integer(is.iter),as.integer(isss),paramfile,as.integer(use_ms_tf_scaling))
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
  list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,sigma = est.cov,ess = ess,pss=pss,params=params,wvec=wvec)		
}


ismoments_parallel = function(site,sigma,mu,is.iter,isss) #numnit=5000 nisamp=20
{
  cd = paste("cd folder", site, ";", sep="")
  folder = paste("folder", site, sep="")
  paramfile <- paste(folder, "/paramfiles_gep/paramfile_new", site, sep="")
  target <- paste(pop, "/gtree_file", site, sep="")
  
  nparams = length(mu) 
  trans.sigma = sigma
  diag(trans.sigma) = 2*diag(trans.sigma)
  params = rmnorm(is.iter, mean = mu, varcov = trans.sigma) #simulate from widened cavity distribution (acts as a prior)
  #use a widened distribution to improve stability
  #correct it below in simwt
  simwt = dmnorm(params,mean=mu,varcov=sigma)/dmnorm(params,mean=mu,varcov=trans.sigma)
  write(t(params), paramfile,ncol=nparams)
  paramfile <- paste("paramfiles_gep/paramfile_new", ilocus, sep="")
  
  
  use_ms_tf_scaling=T #see runnit.r for explanation
  cm1 = paste(cd, "./is_moments_new",target,as.integer(is.iter),as.integer(isss),paramfile,as.integer(use_ms_tf_scaling))
  wvec = system(cm1,intern=T)
  wvec = as.numeric(wvec)
  wvec = simwt*wvec #apply correction
  ess = sum(wvec)^2/sum(wvec^2) #the effective sample size 
  pss = ess/is.iter
  
  if(length(params[,1])!=length(wvec)){
    #print(paste("length params ",length(params[,1])))
    #print(paste("length wvec ",length(wvec)))
    stop("error in length of wvec")
  }
  
  if(ess <= nparams + 2){
    #print(paste("ess problem for",target,": returning prior"))
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
  list(sigmainv.mu=est.sigmainv.mu,sigmainv = est.precmat,mu = est.mu,sigma = est.cov,ess = ess,pss=pss,params=params,wvec=wvec)		
}

loglike <- function(site, params, isss=5000, pop="growing"){
  # Name of file containing phi (x)
  paramfile <- paste("paramfiles/paramfile_new", site, sep="")
  write(params, paramfile, ncol=length(params)) 
  # Create and store file names for target data points of chunk
  target <- paste(pop, "/gtree_file", site, sep="")
  # Now compute the product of all the likelihoods for the chunk
  ll <- log(as.numeric(system(paste("./is_moments_new",target,as.integer(1),as.integer(isss),paramfile, T), intern=T)))
  return(ll)
}

rwmh <- function(start, niter, logtarget, Sigma=NULL){
  # Starting Point and Its Density
  z <- start
  pz <- logtarget(z)
  # Get dimensionality of z when column or row vector
  if (is.null(nrow(z))){
    d <- length(z)
  } else {
    d <- nrow(z)
  }
  # Store Samples
  samples <- matrix(0, nrow=niter, ncol=d)
  samples[1, ] <- z
  # Uniform Random Numbers
  log_u <- log(runif(niter))
  # Proposal: Normal Distribution
  if (is.null(Sigma)){
    Sigma <- diag(d)
  }
  normal_shift <- mvrnorm(n=niter, mu=rep(0, d), Sigma=Sigma)
  for (i in 2:niter){
    # Sample Candidates
    candidate <- z + normal_shift[i, ]
    # Candidate Densities
    p_candidate <- logtarget(candidate)
    # Accept/Reject for Chain 1
    if (log_u[i] <= p_candidate - pz){
      z <- candidate
      pz <- p_candidate
    }
    # Store samples
    samples[i, ] <- z
  }
  return(samples)
}

center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

ismoments2 <- function(il, sigma, mu, Q, r, is.iter, isss){
  # Define tilted distribution (il is the site)
  logtilted <- function(x){
    as.double(-0.5*t(x)%*%(Q%*% x) + t(r) %*% x) + loglike(il, x, pop=pop, isss=isss)
  }
  samples <- rwmh(start = mu, niter = (is.iter+1), logtarget = logtilted)[2:(is.iter+1), ]
  # find qr decomposition
  R_term <- qr.R(qr(center_colmeans(samples)))
  if (any(diag(R_term) == 0)) return(list(sigmainv.mu=r, sigmainv=Q))
  Q_tilted <- chol2inv(R_term) * (is.iter - 3 - 2)
  r_tilted <- Q_tilted %*% apply(samples, 2, mean)
  return(list(sigmainv.mu=r_tilted, sigmainv=Q_tilted))
}

to_exp <- function(mu, Sigma){
  # Change to log
  mu_exp    <- c(exp(mu + 0.5*diag(Sigma)))
  sigma_exp <- diag(mu_exp) %*% (exp(Sigma) - 1) %*% diag(mu_exp)
  return(list(mu=matrix(mu_exp), Sigma=sigma_exp))
}

to_exp2 <- function(r, Q){
  mS <- muSigma(r, Q)
  mu_exp    <- c(exp(mS$mu + 0.5*diag(mS$Sigma)))
  Sigma_exp <- diag(mu_exp) %*% (exp(mS$Sigma) - 1) %*% diag(mu_exp)
  return(list(mu=mu_exp, Sigma=Sigma_exp))
}

####################################### SETTINGS
Niter = 10
Npara = 3  
Nloci = 100
pop = "growing"
# nsim_loci = 100
# if(Nloci > nsim_loci) stop("Nloci > nsim_loci")
nisamp = 20
numnit = 5000
Smooth = 20 
if(Smooth < 1)Smooth = 1
#is.iter = 50
report_counter = 1
mainiter = 1
out_file_name = paste(pop, "_mark.txt", sep="")
################################################# PARAMETERS
# Prior
Prior.sigmainv = diag(Npara)*(1/6.0) 
Prior.sigmainv.mu = rep(0,Npara)
# Local
SigmaInv.mu.sum = rep(0,Npara) 
SigmaInv.sum = diag(Npara)*0 
output = matrix(nrow=Niter,ncol=(Npara+Npara*Npara))
#################################################### MAIN LOOP
abctable = NULL
for(iter in 1:Niter){
  # Cavity
  SigmaInv.mu.cav = (Nloci-1)/Nloci*SigmaInv.mu.sum + Prior.sigmainv.mu
  SigmaInv.cav = (Nloci-1)/Nloci*SigmaInv.sum + Prior.sigmainv
  sigma.cav = chol2inv(chol(SigmaInv.cav))
  mu.cav = sigma.cav%*%SigmaInv.mu.cav 
  # Local
  lambda1.sum = 0
  lambda2.sum = 0
  ess.av = 0
  # function to parallelize
  wrapper <- function(i){
    simres <- ismoments_parallel(il, sigma.cav, as.numeric(mu.cav), numnit, nisamp)
    lambda1.sum.delta <- simres$sigmainv.mu - SigmaInv.mu.cav
    lambda2.sum.delta <- simres$sigmainv - SigmaInv.cav
    return(list(l1=lambda1.sum.delta, l2=lambda2.sum.delta))
  }
  
  out_parallel <- mclapply(1:Nloci, wrapper, mc.cores = 6)
  for (il in 1:Nloci){
    lambda1.sum <- lambda1.sum + out_parallel[[il]]$l1
    lambda2.sum <- lambda2.sum + out_parallel[[il]]$l2
  }
  # Update
  #ess.av = ess.av/Nloci
  SigmaInv.mu.sum = (Smooth - 1)/Smooth*SigmaInv.mu.sum + lambda1.sum/Smooth
  SigmaInv.sum = (Smooth - 1)/Smooth*SigmaInv.sum + lambda2.sum/Smooth
  # Reporting
  sigmat = chol2inv(chol(SigmaInv.sum +  Prior.sigmainv))
  muvec = as.numeric(sigmat%*%(SigmaInv.mu.sum + Prior.sigmainv.mu))
  #sdev = sqrt(diag(sigmat)) #standard deviation
  #smat = diag(1/sdev)
  #corrmat = smat%*%sigmat%*%smat #gets correlation matrix rather than covariance (easier to visualise)
  #index1 = upper.tri(corrmat)
  expparam <- to_exp(muvec, sigmat)
  cat("ExpMu: ", expparam$mu, "ExpSigma: ", diag(expparam$Sigma), "\n")
  report1 = c(muvec,c(sigmat)) #just reports this as a row of output (note the correlation matrix is flattened)
  
  if(report_counter==1){
    write(c(iter,report1),file=out_file_name,ncol=length(report1)+1,append=F)
  }else{
    write(c(iter,report1),file=out_file_name,ncol=length(report1)+1,append=T)
  }
  report_counter = report_counter + 1
  output[iter,] = report1
}
final.out = apply(output[c((ceiling(0.8*Niter)):Niter),],2,mean)

