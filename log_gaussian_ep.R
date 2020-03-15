library(emdbook)   # For dmvnorm
library(MASS)      # For rmvnorm
library(parallel)  # For mclapply

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
sample_moments_tilted <- function(site_index, local_likelihood, r_c, Q_c, nsamples){
  # Create wrapper for tilted distribution. Will be used to sample from it
  # TODO: Check dimensions of input x and output. I use as.double() cause its a vector
  log_target <- function(x){
    log(local_likelihood(x)) + as.double(-0.5 * x %*% (Q_c %*% x) + t(r_c) %*% x) 
  }
  # Grab mu of the cavity distribution as a starting point
  mu_c <- solve(Q_c, r_c)
  # Sample from the targe distribution. First, find the mode and start sampling from there
  sol <- optim(c(mu_c), function(x) -log_target(x), method="BFGS", hessian=TRUE)
  inv_hessian <- chol2inv(chol(sol$hessian))
  samples <- rwmh(start=sol$par, niter=nsamples, logtarget=log_target, Sigma=inv_hessian)
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
  # By sampling tilted, find new global natural parameters
  rq_new <- sample_moments_tilted(site_index, local_likelihood, r_c, Q_c, nsamples)
  # Find new local natural parameters by subtracting natural params of cavity distribution
  r_i_new <- rq_new$r - ep$r + matrix(ep$r_list[,,site_index])
  q_i_new <- rq_new$Q - ep$Q + ep$Q_list[,,site_index]
    # Update EP local parameters
    # Return the new local natural parameter so we can do updates in a parallel fashion
    return(list(r_i=r_i_new, 
                q_i=q_i_new))
}

### Main function for EP Algorithm
ep_algorithm <- function(nsites, nparams, nsamples, npasses, mu_prior, Sigma_prior, local_likelihood, alpha=1.0){
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
    #updates <- mclapply(1:nsites, site_update_wrapper, mc.cores = 4)
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
  }
  # TODO: Should I move this inside the main loop for printing reasons?
  global_moments <- natural2moment(ep$r, ep$Q)
  ep$mu    <- global_moments$mu
  ep$Sigma <- global_moments$Sigma
  return(ep)
}


### EXAMPLE (LIKELIHOOD IS NORMAL)
# Define a likelihood being N((1,1), I)
local_likelihood <- function(x){
  return(dmvnorm(x, mu=c(1, 1), Sigma=diag(2)))
}

# Settings
nsites <- 20
nparams <- 2
nsamples <- 1000
npasses <- 5
mu_prior <- matrix(c(0,0))
Sigma_prior <- diag(2)
alpha <- 1.0

# Run Algorithm
a <- ep_algorithm(nsites, nparams, nsamples, npasses, 
                  mu_prior, Sigma_prior, local_likelihood, alpha)


