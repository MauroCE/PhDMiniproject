library(MASS)
library(pracma)
library(parallel)
library(mnormt)


center_colmeans <- function(x) {
  xcenter = colMeans(x)
  x - rep(xcenter, rep.int(nrow(x), ncol(x)))
}

### RANDOM WALK METROPOLIS HASTINGS
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



ismoments = function(site,sigma,mu,is.iter,isss) 
{
  # command
  cd = paste("cd folder", site, ";", sep="")
  folder = paste("folder", site, sep="")
  paramfile <- paste(folder, "/paramfiles_gep/paramfile_new", site, sep="")
  target <- paste(pop, "/gtree_file", site, sep="")
  # mark stuff
  nparams = length(mu) 
  trans.sigma = sigma
  diag(trans.sigma) = 2*diag(trans.sigma)
  params = rmnorm(is.iter, mean = mu, varcov = trans.sigma) #simulate from widened cavity distribution (acts as a prior)
  #use a widened distribution to improve stability
  #correct it below in simwt
  simwt = dmnorm(params,mean=mu,varcov=sigma)/dmnorm(params,mean=mu,varcov=trans.sigma)
  #paramfile <- paste("paramfiles/paramfile_new", ilocus, sep="")
  #write(t(params), paramfile,ncol=nparams)
  write(params, paramfile, ncol=length(params)) 
  paramfile <- paste("paramfiles_gep/paramfile_new", site, sep="")
  
  use_ms_tf_scaling=T #see runnit.r for explanation
  cm1 = paste(cd, "./is_moments_new",target,as.integer(is.iter),as.integer(isss),paramfile,as.integer(use_ms_tf_scaling))
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



is_positive_definite <- function(Sigma, tol=1e-12){
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  bool <- !all(ev >= -tol*abs(ev[1L]))
  ifelse(bool, F, T)
}

### LOG LIKELIHOOD FUNCTION
loglike <- function(site, params, isss=5000, pop="stable"){
  # cd command
  cd = paste("cd folder", site, ";", sep="")
  folder = paste("folder", site, sep="")
  # Name of file containing phi (x)
  paramfile <- paste(folder, "/paramfiles_gep/paramfile_new", site, sep="")
  write(params, paramfile, ncol=length(params)) 
  paramfile <- paste("paramfiles_gep/paramfile_new", site, sep="")
  # Create and store file names for target data points of chunk
  target <- paste(pop, "/gtree_file", site, sep="")
  # Now compute the product of all the likelihoods for the chunk
  ll <- log(as.numeric(system(paste(cd ,"./is_moments_new",
                                    target,as.integer(1),as.integer(isss),paramfile, T), intern=T)))
  return(ll)
}

to_exp2 <- function(r, Q){
  S <- chol2inv(chol(Q))
  m <- S %*% r
  m_exp    <- c(exp(m + 0.5*diag(S)))
  S_exp <- diag(m_exp) %*% (exp(S) - 1) %*% diag(m_exp)
  return(list(mu=m_exp, Sigma=S_exp))
}

# Function for one sweep
sweep_loci <- function(site, r_global, Q_global, r_sites, Q_sites, nsamples){
  # Set here the deltas in case of non positive definiteness
  DeltaQ <- matrix(0.0, nrow=3, ncol=3)
  Deltar <- matrix(0.0, nrow=3, ncol=1)
  # Form cavity (r, Q) and (m, S)
  r_cavity <- r_global - matrix(r_sites[,,site])
  Q_cavity <- Q_global - Q_sites[,,site]
  S_cavity <- chol2inv(chol(Q_cavity))
  m_cavity <- S_cavity %*% r_cavity
  # S_cavity positive definite?
  if (!is_positive_definite(S_cavity, tol = 1e-9)) {
    return(list(Deltar=Deltar, DeltaQ=DeltaQ))
  }
  # Sample tilted distribution
  logtilted <- function(x){
    -0.5*t(x) %*% (Q_cavity %*% x) + t(r_cavity) %*% x + loglike(site, x, pop=pop)
  }
  # discard first half of the samples
  #nsamples_burnin <- floor(nsamples / 2)
  #samples <- rwmh(start=m_cavity, niter=(nsamples+1), logtarget=logtilted, Sigma=S_cavity)[2:(nsamples+1), ]
  #samples <- rwmh(start=m_cavity, niter=(nsamples+1), logtarget=logtilted, Sigma=S_cavity)[nsamples_burnin:(nsamples+1), ]
  is_out <- ismoments(site,S_cavity, as.numeric(m_cavity),nsamples,5000) 
  # QR decomp + ESS + positive definiteness check
  #Q_tilted <- chol2inv(qr.R(qr(center_colmeans(samples)))) * (nsamples - 3 - 2)
  #r_tilted <- Q_tilted %*% matrix(apply(samples, 2, mean))
  Q_tilted <- is_out$sigmainv   
  r_tilted <- matrix(is_out$sigmainv.mu)
  if (is_positive_definite(Q_tilted)){
    DeltaQ <- Q_tilted - Q_cavity - Q_sites[,,site]
    Deltar <- r_tilted - r_cavity - matrix(r_sites[,,site])
  }
  return(list(Deltar=Deltar, DeltaQ=DeltaQ))
}

####################################### SETTINGS
niter       <- 5
nsamples    <- 100
nfiles      <- 70
ncores      <- 7
mu_prior    <- matrix(rep(0, 3))
Sigma_prior <- diag(3)
pop         <- "growing"
alpha       <- 0.05    # damping factor
setwd("/home/mauro/Documents/University/PhDMiniproject/parallelizing")

################################################# RUN
### INITIALIZATION
# Transform to Q and r parameters
Q_prior <- chol2inv(chol(Sigma_prior))
r_prior <- Q_prior %*% mu_prior
# Store local parameters for each loci (100)
r_sites  <- array(0.0, dim=c(3, 1, nfiles))
Q_sites  <- array(0.0, dim=c(3, 3, nfiles))
# Form global approximation parameter (without the prior, it will be added at the end!)
r_global <- r_prior #apply(r_sites, c(1, 2), sum) + r_prior
Q_global <- Q_prior #apply(Q_sites, c(1, 2), sum) + Q_prior
# MAIN LOOP
for (iter in 1:niter){
  cat("### Iteration: ", iter, "\n")
  sweep_out <- mclapply(1:nfiles, function(i) sweep_loci(i, r_global, Q_global, r_sites, Q_sites, nsamples), mc.cores=ncores)
  # Update (r, Q) global
  Q_global_tentative <- Q_global
  r_global_tentative <- r_global
  for (site in 1:nfiles){
    Q_global_tentative <- Q_global_tentative + alpha*sweep_out[[site]]$DeltaQ
    r_global_tentative <- r_global_tentative + alpha*sweep_out[[site]]$Deltar
  }
  # Check positive definiteness
  if (!is_positive_definite(Q_global_tentative)){
    cat("Q global not positive definite. Skipping update\n")
    #alpha <- max(alpha*0.9, 0.01)
    next
  }
  Q_global <- Q_global_tentative
  r_global <- r_global_tentative
  for (site in 1:nfiles){
    Q_sites[,,site] <- Q_sites[,,site] + alpha*sweep_out[[site]]$DeltaQ
    r_sites[,,site] <- r_sites[,,site] + alpha*sweep_out[[site]]$Deltar
  }
  alpha <- max(alpha*0.9, 0.01)
  # Printing
  muSigmaExp <- to_exp2(r_global, Q_global)
  cat("Iter: ", iter, "muExp: ", muSigmaExp$mu, "Sigma: ", diag(muSigmaExp$Sigma),"\n")
}



