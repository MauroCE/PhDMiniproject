# GAUSSIAN AVERAGED EXPECTATION PROPAGATION
library(MASS)
#install.packages(list.files(pattern="tar.gz"), repos = NULL)
library(mnormt)
library(pracma)

setwd("/home/mauro/Documents/University/PhDMiniproject/parallelizing2/")


ismoments = function(target,cd,paramfile, paramfile2,sigma,mu,is.iter,isss) 
{
  nparams = length(mu) 
  trans.sigma = sigma
  diag(trans.sigma) = 2*diag(trans.sigma)
  params = rmnorm(is.iter, mean = mu, varcov = trans.sigma) #simulate from widened cavity distribution (acts as a prior)
  #use a widened distribution to improve stability
  #correct it below in simwt
  simwt = dmnorm(params,mean=mu,varcov=sigma)/dmnorm(params,mean=mu,varcov=trans.sigma)
  write(t(params), paramfile,ncol=nparams)

  use_ms_tf_scaling=T #see runnit.r for explanation
  cm1 = paste(cd, "./is_moments_new",target,as.integer(is.iter),as.integer(isss),paramfile2,as.integer(use_ms_tf_scaling))
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

ismoments2 <- function(target,cd,paramfile, paramfile2,sigma,Q, r,mu,is.iter,isss, burnin){
  # Define log tilted
  logtilted <- function(x){
    write(x, paramfile, ncol=length(x))
    cm1 = paste(cd, "./is_moments_new",target,as.integer(1),as.integer(isss),paramfile2,as.integer(T))
    as.double(-0.5*t(x) %*% (Q %*% x) + t(r) %*% x) + log(as.numeric(system(cm1, intern=T)))
  }
  # Sample
  samples <- rwmh(start=matrix(mu), niter=(is.iter+1), logtarget=logtilted, burnin = burnin)[2:(is.iter+1), ]
  R <- qr.R(qr(center_colmeans(samples)))
  if (any(diag(R) == 0)) return(list(sigmainv.mu=r, sigmainv=Q))
  Q_tilted <- chol2inv(R) * (is.iter - 3 - 2)
  r_tilted <- Q_tilted %*% apply(samples, 2, mean)
  return(list(sigmainv.mu=r_tilted, sigmainv=Q_tilted))
}


center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
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

rwmh <- function(start, niter, logtarget, burnin=0, Sigma=NULL){
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
  log_u <- log(runif(niter+burnin))
  # Proposal: Normal Distribution
  if (is.null(Sigma)){
    Sigma <- diag(d)
  }
  normal_shift <- mvrnorm(n=(niter+burnin), mu=rep(0, d), Sigma=Sigma)
  for (i in 2:(niter+burnin)){
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
    if (i >= (1+burnin)){
      samples[i-burnin, ] <- z
    }
  }
  return(samples)
}

####################################### SETTINGS
Niter = 10#100
Npara = 3  
Nloci = 60
pop = "growing"
# nsim_loci = 100
# if(Nloci > nsim_loci) stop("Nloci > nsim_loci")
nisamp = 500  #20
n_burn = 80
numnit = 200    #5000
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
  # wrapper function
  wrapper <- function(i){
    # file stuff
    folder <- paste("folder", i, sep="")
    cd <- paste("cd ~/Documents/University/PhDMiniproject/parallelizing2/", folder, " ;", sep="")
    target    <- paste(pop, "/gtree_file", i, sep="")
    paramfile <- paste(folder, "/paramfiles/paramfile_new", i, sep="")
    paramfile2 <- paste("paramfiles/paramfile_new", i, sep="")
    # run mark stuff
    simres <- ismoments2(target,cd,paramfile,paramfile2,sigma.cav,SigmaInv.cav, SigmaInv.mu.cav,as.numeric(mu.cav),numnit,nisamp, n_burn)
    return(list(d1=(simres$sigmainv.mu - SigmaInv.mu.cav), d2=(simres$sigmainv - SigmaInv.cav)))
  }
  out_parallel <- parallel::mclapply(1:Nloci, wrapper, mc.cores = 6)
  for (il in 1:Nloci){
    lambda1.sum <- lambda1.sum + out_parallel[[il]]$d1
    lambda2.sum <- lambda2.sum + out_parallel[[il]]$d2
  }
  # Sweep through all loci
  # for(il in 1:Nloci){
  #   # Sample tilted
  #   ilocus = il
  #   target = paste(pop, "/gtree_file",ilocus,sep="")
  #   simres = ismoments(target,sigma.cav,as.numeric(mu.cav),numnit,nisamp) 
  #   # Get tilted parameters
  #   lambda1.sum = lambda1.sum + simres$sigmainv.mu - SigmaInv.mu.cav
  #   lambda2.sum = lambda2.sum + simres$sigmainv - SigmaInv.cav
  #   ess.av = ess.av + simres$ess
  # }
  # Update
  #ess.av = ess.av/Nloci
  SigmaInv.mu.sum = (Smooth - 1)/Smooth*SigmaInv.mu.sum + lambda1.sum/Smooth
  SigmaInv.sum = (Smooth - 1)/Smooth*SigmaInv.sum + lambda2.sum/Smooth
  # Reporting
  sigmat = try(chol2inv(chol(SigmaInv.sum +  Prior.sigmainv)), silent=TRUE)
  if (class(sigmat) != "matrix") sigmat <- pinv(SigmaInv.sum +  Prior.sigmainv)
  muvec = as.numeric(sigmat%*%(SigmaInv.mu.sum + Prior.sigmainv.mu))
  #sdev = sqrt(diag(sigmat)) #standard deviation
  #smat = diag(1/sdev)
  #corrmat = smat%*%sigmat%*%smat #gets correlation matrix rather than covariance (easier to visualise)
  #index1 = upper.tri(corrmat)
  report1 = c(muvec,c(sigmat)) #just reports this as a row of output (note the correlation matrix is flattened)
  
  # if(report_counter==1){
  #   write(c(iter,ess.av,ess.av/numnit,report1),file=out_file_name,ncol=length(report1)+4,append=F)
  # }else{
  #   write(c(iter,ess.av,ess.av/numnit,report1),file=out_file_name,ncol=length(report1)+4,append=T)
  # }
  # report_counter = report_counter + 1
  output[iter,] = report1
}
final.out = apply(output[c((ceiling(0.8*Niter)):Niter),],2,mean)
to_exp(matrix(final.out[1:3]), matrix(final.out[4:12], 3, 3))

expfunc <- function(m){
  return(to_exp(matrix(m[1:3]), matrix(m[4:12], 3, 3)))
}
