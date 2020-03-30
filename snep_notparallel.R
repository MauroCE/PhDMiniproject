library(ggplot2)
library(tidyr)
### MAIN ALGORITHM
# t: Total number of iterations
snep_async <- function(nworkers, nouter, nsync, mu_prior, Sigma_prior, mu_local, Sigma_local, betas, like,
                 epsilon, t=100, convergence=FALSE, tol=1e-6, maxiter=100){
  # Dimensionality of parameter 
  d <- nrow(Sigma_prior)
  # Steps: 2 to 5
  theta_0  <- to_natural(mu_prior, Sigma_prior)    # Natural parameter of prior distribution
  lambda_i <- to_natural(mu_local, Sigma_local)    # Natural parameter of local approximation
  gamma_i  <- to_mean(mu_local, Sigma_local)       # Mean parameter of local approximation
  theta_c  <- theta_0 + (nworkers - 1)*lambda_i      # Natural parameter of cavity distribution
  theta_p  <- theta_c + lambda_i                   # Auxiliary Parameter: Natural parameter of local global approximation
  # Sample start from a MVN with mean, vcov given by auxiliary parameter
  Sigma_p <- solve(-2*matrix(theta_p[(d+1):(d^2+d), ], nrow=d))
  mu_p    <- Sigma_p %*% matrix(theta_p[1:d, ], ncol=1)
  starts <- mvrnorm(n=nworkers, mu = mu_p, Sigma = Sigma_p)
  # Initialize the server (call it theta_posterior cause it's the on ly field that matters)
  theta_posterior <-  theta_0 + nworkers*lambda_i    # This should be equal to theta_p actually
  # Initialize workers
  workers <- list()
  for (i in 1:nworkers){
    # Store the starting x_i in the history of samples
    x_history      <- matrix(0.0, nrow=max(t+1, maxiter+1), ncol=d)
    x_history[1, ] <- starts[i, ]
    workers[[i]] <- list(x          = matrix(starts[i, ]),
                         gamma      = gamma_i,
                         lambda     = lambda_i,
                         lambda_old = lambda_i,
                         theta_p    = theta_p,
                         beta       = betas[i],
                         theta_c    = theta_c,
                         delta      = 0.0,              # Initialize the delta at zero
                         i          = i,                # Index, used for likelihood
                         x_history  = x_history
                         )                
  }
  # Run workers SYNCHRONOUSLY
  if (!convergence){   # Based on provided number of iterations t
    for (iter in 1:t){
      workers <- run_all_workers_once(workers, iter)
      if (iter %% nouter == 0) {workers <- update_auxiliary(workers)}
      if (iter %% nsync  == 0) {
        comm <- communicate(workers, theta_posterior)
        workers <- comm$w
        theta_posterior <- comm$s
      }
    } 
  } else {
    iter <- 1   # Keep track of iteration number anyways
    distance <- Inf
    while (distance > tol){
      workers <- run_all_workers_once(workers, iter)
      if (iter %% nouter == 0) {workers <- update_auxiliary(workers)}
      if (iter %% nsync  == 0) {
        comm <- communicate(workers, theta_posterior)
        # Update distance, to check convergence
        distance <- sqrt(sum((theta_posterior - comm$s)^2))
        workers <- comm$w
        theta_posterior <- comm$s
      }
      iter <- iter + 1     # Update counter
      if (iter > maxiter){
        cat("Max Iteration Reached. Distance is: ", distance)
        break
      }
    }
  }
  return(list(w=workers, s=theta_posterior))
} 


# RUN ALL WORKERS ONCE
run_all_workers_once <- function(workers, iter){
  for (j in 1:length(workers)){
    # Find mu and Sigma of cavity for every worker (used by likelihood)
    natural_c <- natural2mean(workers[[j]]$theta_c, d, returnMuSigma = TRUE)
    # Skip update if Sigma of cavity distribution is not positive definite
    if (!is_positive_definite(natural_c$Sigma)){
      cat("Skipping worker", j, "at iteration", iter, "as Sigma of cavity distribution is not positive definite.\n")
      next
    }
    # Update state x_i
    log_tilted <- function(x){
      first  <- as.double(t(workers[[j]]$theta_p - (workers[[j]]$lambda / workers[[j]]$beta)) %*% s(x))
      second <- log(like(workers[[j]]$i, c(natural_c$mu), natural_c$Sigma)) / workers[[j]]$beta
      return(first + second)
    }
    # TODO: Which scale to use here?
    x_i_new <- matrix(rwmh(start=workers[[j]]$x, niter=2, logtarget=log_tilted)[2, ]) 
    workers[[j]]$x                   <- x_i_new
    workers[[j]]$x_history[iter+1, ] <- x_i_new
    # Update gamma_i and lambda_i (mean and natural param of local approximation respectively)
    workers[[j]]$gamma  <- workers[[j]]$gamma + epsilon*(s(workers[[j]]$x) - natural2mean(workers[[j]]$theta_c + workers[[j]]$lambda, d))
    workers[[j]]$lambda <- mean2natural(workers[[j]]$gamma, d)
  }
  return(workers)
}

# UPDATE AUXILIARY PARAMETER THETA_I' FOR ALL WORKERS AT OUTER ITERATION
update_auxiliary <- function(workers){
  for (j in 1:length(workers)){
    workers[[j]]$theta_p <- workers[[j]]$theta_c + workers[[j]]$lambda
  }
  return(workers)
}

# COMMUNICATE
communicate <- function(workers, theta_posterior){
  for (j in 1:length(workers)){
    # Update Delta_i
    workers[[j]]$delta      <- workers[[j]]$lambda - workers[[j]]$lambda_old   # Send this to posterior
    workers[[j]]$lambda_old <- workers[[j]]$lambda                  # Keep track of last natural parameter sent
    # Send to posterior, to update theta_posterior
    theta_posterior  <- theta_posterior + workers[[j]]$delta
    # Send new theta_posterior to worker, so it can update natural parameter of cavity distribution
    workers[[j]]$theta_c    <- theta_posterior - workers[[j]]$lambda_old
  }
  return(list(w=workers, s=theta_posterior))
}

### HELPER FUNCTIONS
to_natural <- function(mu, Sigma){    # Transforms (mu, Sigma) into a single natural parameter 
  Q <- solve(Sigma)         # Could use chol2inv(chol()) but requires positive definiteness, and it breaks
  return(rbind(Q %*% mu, matrix(-0.5*Q, ncol=1)))   # TODO: Technically t(Q) but Q should be symmetric anyways?
}
to_mean <- function(mu, Sigma){     # Transforms (mu, Sigma) into the mean parameters for a Gaussian distribution
  return(rbind(mu, matrix(Sigma + tcrossprod(mu), ncol=1))) # TODO: Technically t(Sigma + tcrossprod(mu)) but symmetric
}
natural2mean <- function(theta, d, returnMuSigma=FALSE){   # Transforms a complete natural parameter to mean parameter for Gaussian Approx
  Sigma <- -0.5*solve(matrix(theta[(d+1):(d+d^2), ], nrow=d))  # as (-2)^{-1} = -0.5
  mu <- Sigma %*% matrix(theta[1:d, ])
  if (!returnMuSigma){
    return(to_mean(mu, Sigma))
  } else {
    return(list(mu=mu, Sigma=Sigma)) 
  }
}
mean2natural <- function(nu, d){  # Given a (d + d^2) dimensional mean param, it returns its correponding natural param
  mu    <- matrix(nu[1:d, ])
  Sigma <- matrix(nu[(d+1):(d+d^2), ], nrow=d) - tcrossprod(mu)
  return(to_natural(mu, Sigma))
}
s <- function(x) rbind(x, matrix(tcrossprod(x), ncol=1)) # Calculate sufficient statistics
is_positive_definite <- function(Sigma, tol=1e-10){
  # Positive definite matrices need to be symmetric
  if (!isSymmetric.matrix(Sigma)){
    return(FALSE) # If it's not symmetric, cant be positive definite
  }
  # Compute eigenvalues
  eigenvalues <- eigen(Sigma, only.values = TRUE)$values
  # If smaller then threshold, set to zero
  eigenvalues[abs(eigenvalues) < tol] <- 0.0
  # If any is <= 0 then it's not positive definite
  return(ifelse(any(eigenvalues <= 0), FALSE, TRUE))
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
    Sigma <- 5*diag(d)
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


### LIKELIHOOD
like <- function(i, mu, sigma, isss=5000){    # w is a worker
  # Construct file name containing chunk of data
  target <- paste("genetree/gtree_file", i, sep="")
  # Construct name of the file containing the sample from the cavity distribution
  paramfile <- paste("paramfiles/paramfile_new", i, sep="")
  # mu    : cavity distribution mean - prior mean for log(theta) and log(rho)
  # sigma : cavity distribution vcov - prior vcov for log(theta) and log(rho)
  nparams <- length(mu) # = 3 in this case with log(theta), log(R), log(tf), 
  # simulate from widened cavity distribution (acts as a prior). Could use lgarch::rmnorm
  trans.sigma <-  sigma
  diag(trans.sigma) <-  2*diag(trans.sigma)
  params = MASS::mvrnorm(n=1, mu=mu, Sigma=trans.sigma)
  # use a widened distribution to improve stability
  # correct it below in simwt
  simwt = dmvnorm(params, mu=mu, Sigma=sigma) / dmvnorm(params, mu=mu, Sigma=trans.sigma)
  write(t(params),paramfile,ncol=nparams)
  use_ms_tf_scaling=T #see runnit.r for explanation
  cm1 = paste("./is_moments_new",target,as.integer(1),as.integer(isss),paramfile,as.integer(use_ms_tf_scaling))
  wvec = system(cm1,intern=T)
  wvec = as.numeric(wvec)
  wvec = simwt*wvec #apply correction
  return(wvec)
}


### RUN EXAMPLE
d           <- 3                           # Dimensionality of parameter space. Here 3 cause (logTheta, logR, logTf) 
nworkers    <- 2
nouter      <- 8
nsync       <- 8
t           <- 200
mu_prior    <- mu_local    <- matrix(rep(0, d))
Sigma_prior <- Sigma_local <- diag(d)
betas       <- rep(1/nworkers, nworkers)   # pSNEP
epsilon     <- 0.02                        # learning rate
convergence <- FALSE                        # Whether to stop algorithm based on number of iterations t or on convergence
tol         <- 1e-6                        # Tolerance for checking convergence
maxiter     <- 1000
out <- snep_async(nworkers=nworkers, 
                  nouter=nouter, 
                  nsync=nsync, 
                  mu_prior=mu_prior, 
                  Sigma_prior=Sigma_prior, 
                  mu_local=mu_local, 
                  Sigma_local=Sigma_local,
                  betas=betas,
                  epsilon=epsilon,
                  like=like,
                  t=t,
                  convergence=convergence,
                  tol=tol,
                  maxiter=maxiter)

exp(natural2mean(out$s, d, TRUE)$mu)
exp(natural2mean(out$s, d, TRUE)$Sigma)


# Prepare data from worker 1
df1       <- data.frame(out$w[[1]]$x_history)
df1$index <- 1:nrow(df1)
df1$worker <- 1
# Prepare data from worker 2
df2 <- data.frame(out$w[[2]]$x_history)
df2$index <- 1:nrow(df2)
df2$worker <- 2
# Melt together
melt <- gather(rbind(df1, df2), "param_dimension", "trace", -index, -worker)
ggplot(data=melt) +
  geom_line(aes(x=index, y=trace)) +
  facet_wrap(worker~param_dimension)
