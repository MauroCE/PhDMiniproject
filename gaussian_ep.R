# GAUSSIAN EXPECTATION PROPAGATION
library(mcmc)
library(emdbook)

###########################
### Normal Distribution ###
###########################
# Multivariate Normal Density function specified
# using the natural parameters r and Q. It is unnormalized
# Since I'm not coding A(r, Q)
## r : Sigma^{-1} * mu
## Q : Sigma^{-1}
normal <- function(x, r, Q){
  return(exp(r %*% x - 0.5 * x %*% (Q %*% x)))
}
# The following is the same as the above function, but on the 
# log scale. This should help computations and sampling
log_normal <- function(x, r, Q){
  return(r %*% x - 0.5 * x %*% (Q %*% x))
}

#############################
### Converting Parameters ###
#############################
# From "moment" parameters (mu, Sigma) to natural parameters
# i.e. (r, Q). Technically would be -0.5*Q, but practically this doesn't matter
# as it can be factored out.
moment2natural <- function(mu, Sigma){
  Q <- solve(Sigma) #chol2inv(chol(Sigma))   #TODO: CHECK THIS OUT
  r <- Q %*% mu
  return(list(r=r, Q=Q))
}
# From natural parameters (r, Q) to "moment" parameters (mu, Sigma).
natural2moment <- function(r, Q){
  Sigma <- solve(Q) #chol2inv(chol(Q))
  mu <- Sigma %*% r
  return(list(mu=mu, Sigma=Sigma))
}


######################
### Initialization ###
######################
# Function that initializes EP
initialize_ep <- function(nsites, mu_prior, Sigma_prior){
  # Find dimensionality of parameters
  nparams <- nrow(Sigma_prior)
  # Store r, Q for the prior
  natural_prior <- moment2natural(mu_prior, Sigma_prior)
  
  ep <- list(
    r=natural_prior$r, # Set global r and global Q to r_prior and q_prior. 
    Q=natural_prior$Q,
    mu=mu_prior,       # Set global mu and Sigma to mu_prior and Sigma_prior
    Sigma=Sigma_prior,
    r_list=array(0, dim=c(nparams, 1, nsites)),  # Store arrays containing r_i and Q_i for every site
    Q_list=array(0, dim=c(nparams, nparams, nsites))
  )
  return(ep)
}

######################
####### Updates ######
######################
# Checks if the matrix is positive definite. No need to check if it is 
# square or symmetric because by construction they will be
is_positive_definite <- function(Sigma, tol=1e-8){
  # Compute eigenvalues
  eigenvalues <- eigen(Sigma, only.values = TRUE)$values
  # If smaller then threshold, set to zero
  eigenvalues[abs(eigenvalues) < tol] <- 0.0
  # If any is <= 0 then it's not positive definite
  return(ifelse(any(eigenvalues <= 0), FALSE, TRUE))
}


# Finds mu and Sigma of the cavity distribution for site i
moments_cavity <- function(ep, site_index){
  # Subtract Q_i from Q_global, and r_i from r_global to find natural parameters
  # of cavity distribution
  Q_c <- ep$Q - ep$Q_list[,,site_index]
  r_c <- ep$r - ep$r_list[,,site_index]
  moment_cavity <- natural2moment(r_c, Q_c)
  return(list(mu_c=moment_cavity$mu, 
              Sigma_c=moment_cavity$Sigma))
}

# Compute samples from the tilted distribution and then computes sample
# moments from them. 
sample_moments_tilted <- function(site_index, local_likelihood, mu_c, Sigma_c, nsamples){
  # Create wrapper for tilted distribution. Will be used to sample from it
  target <- function(x) local_likelihood(x) * dmvnorm(x, c(mu_c), Sigma_c)
  # Sample from the targe distribution. First, find the mode and start sampling from there
  sol <- optim(rep(0, nrow(Sigma_c)), function(x) -log(target(x)), method="BFGS", hessian=TRUE)
  samples <- metrop(target, sol$par, nbatch=nsamples)$batch
  # Compute mean and covariance of the samples, assign them to the new global moment parameters
  mu_global_new    <- matrix(apply(samples, 2, mean))
  Sigma_global_new <- cov(samples)
  # Convert new global moment parameters to natural parameters
  Q_global_new <- chol2inv(chol(Sigma_global_new))
  r_global_new <- Q_global_new %*% mu_global_new
  return(list(Q=Q_global_new, r=r_global_new))
}

site_update <- function(ep, site_index, pass_index, nsamples, local_likelihood){
  # Compute mu, Sigma of cavity distribution
  moment_cavity <- moments_cavity(ep, site_index)
  Sigma_c <- moment_cavity$Sigma
  mu_c    <- moment_cavity$mu
  # Check for positive definiteness of Sigma, otherwise, skip this
  if (is_positive_definite(Sigma_c)){
    # By sampling tilted, find new global natural parameters
    rq_new <- sample_moments_tilted(site_index, local_likelihood, mu_c, Sigma_c, nsamples)
    # Find new local natural parameters by subtracting natural params of cavity distribution
    r_i_new <- rq_new$r - ep$r + matrix(ep$r_list[,,site_index])
    q_i_new <- rq_new$Q - ep$Q + ep$Q_list[,,site_index]
    # Update EP local parameters
    # TODO: Do this after everything happens, to be thredsafe
    #ep$r_list[,,site_index] <- r_i_new
    #ep$Q_list[,,site_index] <- q_i_new
    # Return the new local natural parameter so we can do updates in a parallel fashion
    return(list(r_i=r_i_new, 
                q_i=q_i_new))
  } else {
    cat("Skipping site", site_index, "in pass", pass_index, "as Sigma of Cavity is not positive definite.\n")
    #TODO: Do we keep the old ones, or pass empty ones?
    # 
    # return(list(
    #   r_i = matrix(0, ncol=1, nrow=nrow(Sigma_c)),
    #   q_i = matrix(0, nrow=nrow(Sigma_c), ncol=nrow(Sigma_c))
    # ))
    return(list(
      r_i = matrix(ep$r_list[,,site_index]),
      q_i = ep$Q_list[,,site_index]
    ))
  }
}

ep_algorithm <- function(nsites, nparams, nsamples, npasses, mu_prior, Sigma_prior, local_likelihood, alpha=1.0){
  # Store prior natural parameters, arrays of Q_i and r_i, global Q and r
  ep <- initialize_ep(nsites, mu_prior, Sigma_prior)
  # Each pass computes updates in parallel
  for (pass_index in 1:npasses){
    # Create a wrapper function so that we can use mclapply
    site_update_wrapper <- function(site_index){
      return(site_update(ep, site_index, pass_index, nsamples, local_likelihood))
    }
    # Find mu, Sigma of cavity distribution, construct tilted, sample from it and convert moment
    # parameter to natural parameters.
    #TODO: COULD USE MCLAPPLY BUT FOR NOW LOOP THROUGH
    # updates <- mclapply(1:nsites, site_update_wrapper, mc.cores=4)
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
    # Update the global mu and Sigma
    global_moments <- natural2moment(ep$r, ep$Q)
    ep$mu <- global_moments$mu
    ep$Sigma <- global_moments$Sigma
  }
  return(ep)
}


local_likelihood <- function(x){
  return(dmvnorm(x, mu=c(1,1), Sigma=diag(2)))
}

  
# EXAMPLE
nsites <- 20
nparams <- 2
nsamples <- 100
npasses <- 10
mu_prior <- matrix(c(0,0))
Sigma_prior <- diag(2)
alpha <- 0.01
a <- ep_algorithm(nsites, nparams, nsamples, npasses, mu_prior, Sigma_prior, local_likelihood, alpha=alpha)


