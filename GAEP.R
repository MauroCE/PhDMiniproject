# GAUSSIAN AVERAGED EXPECTATION PROPAGATION
library(MASS)
install.packages(list.files(pattern="tar.gz"), repos = NULL)
library(mnormt)

### LOG LIKELIHOOD FUNCTION
loglike <- function(site, params, isss=5000, pop="growing", simwt=1.0){
  # Name of file containing phi (x)
  paramfile <- paste("paramfiles/paramfile_new", site, sep="")
  write(params, paramfile, ncol=length(params)) 
  # Create and store file names for target data points of chunk
  target <- paste(pop, "/gtree_file", site, sep="")
  # Now compute the product of all the likelihoods for the chunk
  ll <- log(simwt*as.numeric(system(paste("./is_moments_new",target,as.integer(1),as.integer(isss),paramfile, T), intern=T)))
  return(ll)
}

matrix_inverse <- function(matrix){
  trial <- try(chol2inv(chol(matrix)), silent=TRUE)
  if (class(trial) == "matrix") return(trial) 
  else{
    trial <- try(solve(matrix), silent=TRUE)
    if (class(trial) == "matrix") return(trial)
    else pinv(matrix)
  }
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

#ismoments = function(target,sigma,mu,is.iter,isss) 


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
Niter = 2
Npara = 3  
Nloci = 20 #100
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
  # Sweep through all loci
  for(il in 1:Nloci){
    # Sample tilted
    ilocus = il
    target = paste(pop, "/gtree_file",ilocus,sep="")
    #simres = ismoments(target,sigma.cav,as.numeric(mu.cav),numnit,nisamp) 
    trans.sigma.cav <- sigma.cav
    diag(trans.sigma.cav) <- 2*diag(sigma.cav)
    start <- mvrnorm(n=1, mu=as.numeric(mu.cav), Sigma=trans.sigma.cav)
    
    # Form log-tilted distribution
    logtilted <- function(x){
      simwt = dmnorm(c(x),mean=c(mu.cav),varcov=sigma.cav)/dmnorm(c(x),mean=c(mu.cav),varcov=trans.sigma.cav)
      return(as.double(-0.5*t(x) %*% (SigmaInv.cav %*% x) + t(SigmaInv.mu.cav) %*% x) + loglike(il, x, pop=pop, simwt=simwt))
    }
    # Sample from the tilted distribution
    samples <- rwmh(start=start, niter=(nisamp+1), logtarget=logtilted, Sigma=sigma.cav)[2:(nisamp+1), ]
    # Obtain empirical estimates of mu_tilted and Sigma_tilted
    mu_tilted    <- apply(samples, 2, mean)
    Sigma_tilted <- cov(samples)
    # transform them to natural params
    simres.sigmainv    <- matrix_inverse(Sigma_tilted)
    simres.sigmainv.mu <- simres.sigmainv %*% mu_tilted
    # Get tilted parameters
    lambda1.sum = lambda1.sum + simres.sigmainv.mu - SigmaInv.mu.cav
    lambda2.sum = lambda2.sum + simres.sigmainv - SigmaInv.cav
    #ess.av = ess.av + simres$ess
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
  report1 = c(muvec,c(sigmat)) #just reports this as a row of output (note the correlation matrix is flattened)
  
  if(report_counter==1){
    #write(c(iter,ess.av,ess.av/numnit,report1),file="growing_mark.txt",ncol=length(report1)+4,append=F)
    write(c(iter,report1),file="growing_mark_mine.txt",ncol=length(report1)+1,append=F)
  }else{
    write(c(iter,report1),file="growing_mark_mine.txt",ncol=length(report1)+1,append=T)
  }
  report_counter = report_counter + 1
  output[iter,] = report1
}
final.out = apply(output[c((ceiling(0.8*Niter)):Niter),],2,mean)

