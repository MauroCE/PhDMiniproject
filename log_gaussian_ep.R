library(emdbook)   # For dmvnorm
library(MASS)      # For mvrnorm
library(parallel)  # For mclapply
# library(lgarch)    # For rmnorm
library(mnormt)    # For dmnorm, rmnorm

# Checks for positive definiteness
is_positive_definite <- function(Sigma, tol=1e-8){
  # Compute eigenvalues
  eigenvalues <- eigen(Sigma, only.values = TRUE)$values
  # If smaller then threshold, set to zero
  eigenvalues[abs(eigenvalues) < tol] <- 0.0
  # If any is <= 0 then it's not positive definite
  return(ifelse(any(eigenvalues <= 0), FALSE, TRUE))
}

### Random-Walk Metropolis-Hastings (on log scale)
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

### Converting from Natural Parameters to Moment Parameters and vice-versa
moment2natural <- function(mu, Sigma){
  Q <- solve(Sigma)  # chol2inv(chol(Sigma)) is faster but sometimes fails
  r <- Q %*% mu
  return(list(r=r, Q=Q))
}
natural2moment <- function(r, Q){
  Sigma <- solve(Q) # chol2inv(chol(Sigma)) is faster but sometimes fails
  mu <- Sigma %*% r
  return(list(mu=mu, Sigma=Sigma))
}

### Sample Moments from Tilted Distribution
sample_moments_tilted <- function(site_index, local_likelihood, mu_c, Sigma_c, nsamples){
  # Create wrapper for tilted distribution. Will be used to sample from it
  # TODO: Check dimensions of input x and output. I use as.double() cause its a vector
  # Find moment parameters for cavity distribution (will be used by likelihood)
  log_tilted <- function(theta){
    # Evaluate log Cavity Distribution (prior) at theta using r_c and Q_c
    # log_cavity <- as.double(-0.5 * theta %*% (Q_c %*% theta) + t(r_c) %*% theta)
    log_cavity <- dmvnorm(theta, c(mu_c), Sigma_c, log=TRUE)
    # Evaluate log Local Likelihood at the data contained at site `site_index`, given r_c and Q_c
    log_likelihood <- local_likelihood(site_index, sigma=Sigma_c, mu=c(mu_c))
    # Sum them log_cavity and log_likelihood to obtain log_tilted
    # Compute Cavity Distribution Estimate
    return(log_likelihood + log_cavity)
  }
  # Sample from the targe distribution. First, find the mode and start sampling from there
  sol <- optim(c(mu_c), function(x) -log_tilted(x), method="BFGS", hessian=TRUE)
  inv_hessian <- chol2inv(chol(sol$hessian))
  samples <- rwmh(start=sol$par, niter=nsamples, logtarget=log_tilted, Sigma=inv_hessian)
  # Compute mean and covariance of the samples, assign them to the new global moment parameters
  mu_global_new    <- matrix(apply(samples, 2, mean))
  Sigma_global_new <- cov(samples)
  # Convert new global moment parameters to natural parameters
  natural_global_new <- moment2natural(mu_global_new, Sigma_global_new)
  return(list(Q=natural_global_new$Q, r=natural_global_new$r))
}

### Performs Site update (Moment Matching)
site_update <- function(ep, site_index, pass_index, nsamples, local_likelihood){
  # Compute natural parameters of cavity distribution
  Q_c <- ep$Q - ep$Q_list[,,site_index]
  r_c <- ep$r - matrix(ep$r_list[,,site_index])
  # Compute moment parameters to see whether to skip this site update
  moment_cavity <- natural2moment(r_c, Q_c)
  if (is_positive_definite(moment_cavity$Sigma)){
    # By sampling tilted, find new global natural parameters
    rq_new <- sample_moments_tilted(site_index, local_likelihood, moment_cavity$mu, moment_cavity$Sigma, nsamples)
    # Find new local natural parameters by subtracting natural params of cavity distribution
    r_i_new <- rq_new$r - ep$r + matrix(ep$r_list[,,site_index])
    q_i_new <- rq_new$Q - ep$Q + ep$Q_list[,,site_index]
    # Update EP local parameters
    # Return the new local natural parameter so we can do updates in a parallel fashion
    return(list(r_i=r_i_new, 
                q_i=q_i_new))
  } else {
    return(list(r_i=matrix(0, nrow=nparams, ncol=1),
                q_i=matrix(0, nrow=nparams, ncol=nparams)))
  }
}

### Main function for EP Algorithm
ep_algorithm <- function(nsites, nparams, nsamples, npasses, mu_prior, Sigma_prior, local_likelihood, alpha=1.0){
  # Create arrays to store the history of the r and Q
  history_r <- array(0, dim=c(nparams, 1, npasses))
  history_Q <- array(0, dim=c(nparams, nparams, npasses))
  # Store prior natural parameters, arrays of Q_i and r_i, global Q and r
  natural_prior <- moment2natural(mu_prior, Sigma_prior)
  ep <- list(
    r=natural_prior$r, # Set global r and global Q to r_prior and q_prior. 
    Q=natural_prior$Q,
    r_list=array(0, dim=c(nparams, 1, nsites)),  # Arrays with r_i and Q_i for every site
    Q_list=array(0, dim=c(nparams, nparams, nsites))
  )
  # Each pass computes updates in parallel
  for (pass_index in 1:npasses){
    cat("Pass", pass_index,": Global r is\n", ep$r, "\nGlobal Q is\n", ep$Q, "\n")
    # Create a wrapper function so that we can use mclapply
    site_update_wrapper <- function(site_index){
      return(site_update(ep, site_index, pass_index, nsamples, local_likelihood))
    }
    # Find (mu, Sigma) of cavity distribution, construct tilted, sample from it 
    # and convert moment parameters to natural parameters.
    # updates <- mclapply(1:nsites, site_update_wrapper, mc.cores = 4)
    updates <- list()
    for (site in 1:nsites){
      updates[[site]] <- site_update_wrapper(site)
    }
    # Sum up all the updates and update ep local natural parameters
    r_new <- 0.0
    Q_new <- 0.0
    for (site in 1:nsites){
      ep$r_list[,,site] <- updates[[site]]$r_i
      ep$Q_list[,,site] <- updates[[site]]$q_i
      r_new <- r_new + updates[[site]]$r_i
      Q_new <- Q_new + updates[[site]]$q_i
    }
    # Update the global parameters with some damping 
    ep$r <- alpha*r_new + (1 - alpha)*(ep$r)
    ep$Q <- alpha*Q_new + (1 - alpha)*(ep$Q)
    # Add global mu and global Sigma to history to check convergence
    history_r[,,pass_index] <- ep$r
    history_Q[,,pass_index] <- ep$Q
  }
  return(list(ep=ep,
              history_r=history_r,
              history_Q=history_Q))
}


### EXAMPLE (LIKELIHOOD IS NORMAL)
# Define a likelihood being N((1,1), I)
# local_likelihood <- function(x){
#   return(dmvnorm(x, mu=c(1, 1), Sigma=diag(2)))
# }

# Notice that can check if a directory exists using 
# dir.exists(file.path(getwd(), "paramfiles"))
# if this doesnt work, one can replace getwd() with the correct "main" path
# which usually is "/home/mauro/Documents/University/MHAAR/IS_sub3"
# could use something like 
# unlist(mclapply(1:nsites, function(i) local_like(i, diag(3), c(0,0,0), 1)))
local_likelihood <- function(site_index, sigma, mu, isss=1000){
  # Construct file name containing chunk of data
  target <- paste("genetree/gtree_file", site_index, sep="")
  # Construct name of the file containing the sample from the cavity distribution
  paramfile <- paste("paramfiles/paramfile_new", site_index, sep="")
  # mu    : cavity distribution mean - prior mean for log(theta) and log(rho)
  # sigma : cavity distribution vcov - prior vcov for log(theta) and log(rho)
  nparams <- length(mu) # = 3 in this case with log(theta), log(R), log(tf), 
  # simulate from widened cavity distribution (acts as a prior). Could use lgarch::rmnorm
  trans.sigma <-  sigma
  diag(trans.sigma) <-  2*diag(trans.sigma)
  params = mnormt::rmnorm(n=1, mean = mu, varcov = trans.sigma) 
  # use a widened distribution to improve stability
  # correct it below in simwt
  simwt = mnormt::dmnorm(params,mean=mu,varcov=sigma) / mnormt::dmnorm(params,mean=mu,varcov=trans.sigma)
  write(t(params),paramfile,ncol=nparams)
  use_ms_tf_scaling=T #see runnit.r for explanation
  cm1 = paste("./is_moments_new",target,as.integer(1),as.integer(isss),paramfile,as.integer(use_ms_tf_scaling))
  wvec = system(cm1,intern=T)
  wvec = as.numeric(wvec)
  wvec = simwt*wvec #apply correction
  return(wvec)
}


# Settings
nsites <- 5
nsamples <- 500
npasses <- 20
mu_prior <- matrix(c(0,0,0))
Sigma_prior <- diag(3)
nparams <- nrow(Sigma_prior)
# For Likelihood
alpha <- 0.2

# Run Algorithm
a <- ep_algorithm(nsites, nparams, nsamples, npasses, 
                  mu_prior, Sigma_prior, local_likelihood, alpha)


# After its done, to look at the history
history_mu    <- array(0, dim=c(nparams, 1, npasses))
history_Sigma <- array(0, dim=c(nparams, nparams, npasses))
for (i in 1:npasses){
  natural <- natural2moment(matrix(a$history_r[,,i]), a$history_Q[,,i])
  history_mu[,,i] <- natural$mu
  history_Sigma[,,i] <- natural$Sigma
}
