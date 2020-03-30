library(ggplot2)
library(tidyr)
library(mcmc)
library(dplyr)
library(parallel)
library(MASS)
### MAIN ALGORITHM
# t: Total number of iterations
snep <- function(nworkers, nouter, nsync, mu_prior, Sigma_prior, mu_local, Sigma_local, betas, loglike,
                 epsilon, pop="stable", tol=1e-6, maxiter=1000){
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
  # Store history of theta_posterior mean
  n_theta_updates <- maxiter %/% nsync
  theta_posterior_history      <- matrix(0.0, nrow=(n_theta_updates+1), ncol=(d+d^2))
  theta_posterior_history[1, ] <- theta_posterior
  # Divide the 20 sites into the correct number of workers
  chunken_sites <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  chunk_sites   <- chunken_sites(1:20, nworkers)
  # Initialize workers
  workers <- list()
  for (i in 1:nworkers){
    # Store the starting x_i in the history of samples
    x_history      <- matrix(0.0, nrow=(maxiter+1), ncol=d)
    x_history[1, ] <- starts[i, ]
    workers[[i]]   <- list(x          = matrix(starts[i, ]),
                           gamma      = gamma_i,
                           lambda     = lambda_i,
                           lambda_old = lambda_i,
                           theta_p    = theta_p,
                           beta       = betas[i],
                           theta_c    = theta_c,
                           delta      = rep(Inf, d+d^2),              # Initialize the delta at Inf
                           sites      = chunk_sites[[i]],             # Indeces of site chunks
                           x_history  = x_history
    )                
  }
  # Run workers SYNCHRONOUSLY
  iter     <- 1
  distance <- sqrt(sum((workers[[1]]$delta)^2))
  while (distance > tol){
    cat("Iteration: ", iter, " Distance: ", distance, " Mu: ", find_mu_from_natural(theta_posterior, d), "\n")
    # RUN ALL WORKERS ONCE (SAMPLING + UPDATING GAMMA AND LAMBDA)
    for (j in 1:nworkers){
      natural_c <- natural2mean(workers[[j]]$theta_c, d, returnMuSigma = TRUE)  # mu, Sigma of cavity distribution
      if (!is_positive_definite(natural_c$Sigma)){                   # Sigma_cavity must be positive definite
        cat("Skipping worker", j, "at iteration", iter, "as Sigma of cavity distribution is not positive definite.\n")
        cat(natural_c$Sigma, "\n")
        next
      }
      # Update state x_i
      # Calculate theta_i' - lambda_i / beta_i outside the log tilted function and calculate 
      tilted_param <- workers[[j]]$theta_p - (workers[[j]]$lambda / workers[[j]]$beta)
      # Calculate log likelihood at the ith site, divided by beta_i  # TODO: ell(x) and s(x), x not the same??
      # ll <- log(like(workers[[j]]$i, c(natural_c$mu), natural_c$Sigma)) / workers[[j]]$beta
      # Log tilted distribution
      log_tilted <- function(x){
        return(as.double(t(tilted_param) %*% s(x)) + loglike(j, workers[[j]]$sites, c(x), pop=pop) / workers[[j]]$beta)
      }   
      x_i_new <- matrix(rwmh(start=workers[[j]]$x, niter=2, logtarget=log_tilted)[2, ]) #TODO: Scale?
      workers[[j]]$x                   <- x_i_new
      workers[[j]]$x_history[iter+1, ] <- x_i_new
      # Update gamma_i and lambda_i (mean and natural param of local approximation respectively)
      workers[[j]]$gamma  <- workers[[j]]$gamma + epsilon*(s(workers[[j]]$x) - natural2mean(workers[[j]]$theta_c + workers[[j]]$lambda, d))
      workers[[j]]$lambda <- mean2natural(workers[[j]]$gamma, d)
    }
    # UPDATE AUXILIARY PARAMETER  (Parallelized)
    if (iter %% nouter == 0) {
      update_auxiliary <- function(j){
        workers[[j]]$theta_p <- workers[[j]]$theta_c + workers[[j]]$lambda
        return(workers[[j]])
      } 
      workers <- mclapply(1:nworkers, update_auxiliary, mc.cores = 4)
    }
    # COMMUNICATE WITH SERVER AND UPDATE THETA_POSTERIOR
    if (iter %% nsync  == 0){
      for (j in 1:nworkers){
        workers[[j]]$delta          <- workers[[j]]$lambda - workers[[j]]$lambda_old   # Send this to posterior
        workers[[j]]$lambda_old     <- workers[[j]]$lambda                  # Keep track of last natural parameter sent
        # Find new theta_posterior
        theta_posterior             <- theta_posterior + workers[[j]]$delta
        workers[[j]]$theta_c        <- theta_posterior - workers[[j]]$lambda_old
      }
      # Update history of theta posterior
      theta_posterior_history[(iter %/% nsync)+1, ] <- theta_posterior
      # The difference between the "old" theta_posterior and the "new" theta_posterior is the sum of the deltas
      total_delta  <- matrix(0, nrow=(d+d^2), ncol=1)
      for (j in 1:nworkers){
        total_delta <- total_delta + workers[[j]]$delta
      }
      distance <- sqrt(sum(total_delta^2))
    }
    iter <- iter + 1
    # Stop if reached maximum number of iterations
    if (iter >= maxiter){
      cat("Max Iteration Reached. Distance is: ", distance)
      break
    }
  }
  # Trim the sample history for every worker, in case that convergence has been achieved earlier than maxiter
  if (iter < maxiter){
    for(j in 1:nworkers){
      workers[[j]]$x_history <- workers[[j]]$x_history[1:iter, ]  # iter, not iter+1 because we've just set iter += 1
    }
  }
  return(list(workers=workers, server=theta_posterior, theta_posterior_history=theta_posterior_history))
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
find_mu_from_natural <- function(theta, d) solve(matrix(-2*theta[(d+1):(d+d^2), ], nrow=d), matrix(theta[1:d, ]))
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


### LOG LIKELIHOOD
loglike <- function(worker_index, sites, params, isss=5000, pop="stable"){
  # Name of file containing phi (x)
  paramfile <- paste("paramfiles/paramfile_new", worker_index, sep="")
  write(t(params),paramfile,ncol=length(params)) 
  # Create and store file names for target data points of chunk
  chunk_files <- list()
  for (i in seq_along(sites)){
    chunk_files[[i]] <- paste(pop, "/gtree_file", sites[[i]], sep="")
  }
  # Now compute the product of all the likelihoods for the chunk
  ll <- 0.0
  for (target in chunk_files){
    ll <- ll + log(as.numeric(system(paste("./is_moments_new",target,as.integer(1),as.integer(isss),paramfile, T), intern=T)))
  }
  return(ll)
}


### RUN EXAMPLE
d           <- 3                                  # Number of parameters. Here 3 cause (logTheta, logR, logTf) 
nworkers    <- 4                                  # Number of sites
nouter      <- 10                                 # Number of iterations after which we update auxiliary parameter
nsync       <- 10                                 # Number of iterations after which we communicate with the server
mu_prior    <- mu_local    <- matrix(rep(0, d))
Sigma_prior <- Sigma_local <- diag(d)
betas       <- rep(1/nworkers, nworkers)          # pSNEP
epsilon     <- 0.01                               # learning rate
pop         <- "stable"                      # Choose between "stable", "growing", or "contracting"
tol         <- 1e-6                               # Tolerance used to determined whether theta_posterior has converged
maxiter     <- 200                                   
# true values
if (pop == "stable"){
  true_values <- c(10, 1, 1)
} else if (pop == "growing"){
  true_values <- c(20, 20, 0.05)
} else {
  true_values <- c(1, 0.05, 1)
}
out <- snep(nworkers=nworkers, 
            nouter=nouter, 
            nsync=nsync, 
            mu_prior=mu_prior, 
            Sigma_prior=Sigma_prior, 
            mu_local=mu_local, 
            Sigma_local=Sigma_local,
            betas=betas,
            epsilon=epsilon,
            loglike=loglike,
            pop=pop,
            tol=tol,
            maxiter=maxiter
)

exp(natural2mean(out$server, d, TRUE)$mu)
exp(natural2mean(out$server, d, TRUE)$Sigma)


# Prepare data from worker 1
plot_traces <- function(){
  create_df <- function(i) data.frame(out$workers[[i]]$x_history) %>% mutate(index=row_number())
  # create names for workers
  worker_names <- c()
  for (j in 1:nworkers) worker_names[paste(j, sep='')] <- paste("worker", j, sep="")
  lapply(1:nworkers, create_df) %>% 
    bind_rows(.id="worker") %>% 
    gather("param_dimension", "trace", -index, -worker) %>% 
    {ggplot(data=., aes(x=index, y=trace)) + 
        geom_line() + 
        facet_wrap(worker ~ param_dimension, labeller = labeller(worker = worker_names))
    }
}
plot_traces()

# Function to see the path of theta_posterior
theta_posterior_path <- function(true_values){
  n <- nrow(out$theta_posterior_history)
  if (all(out$theta_posterior_history[n, ] == rep(0.0, ))){
    tp_history <- out$theta_posterior_history[1:(n-1), ]
    n          <- n-1
  } else {
    tp_history <- out$theta_posterior_history
  }
  means  <- matrix(0.0, nrow=n, ncol=d)
  sds    <- matrix(0.0, nrow=n, ncol=d)
  for (j in 1:n){
    # means
    mean_params <- natural2mean(matrix(tp_history[j, ]), d, T)
    means[j, ] <- exp(mean_params$mu + 0.5*matrix(diag(mean_params$Sigma)))
    sds[j, ]   <- sqrt(diag(diag(means[j, ]) %*% (exp(mean_params$Sigma) - 1) %*% diag(means[j, ])))
  }
  # Find min and max
  min_means <- means - 1.96 * sds
  max_means <- means + 1.96 * sds
  # Get means df
  means_df <- means %>% data.frame()     %>% mutate(index=row_number()) %>% gather("Dimension", "Path", -index) 
  # get min df
  min_df   <- min_means %>% data.frame() %>% mutate(index=row_number()) %>% gather("Dimension", "Path", -index)
  # get max df
  max_df   <- max_means %>% data.frame() %>% mutate(index=row_number()) %>% gather("Dimension", "Path", -index)
  # Final df
  final_df <- mutate(means_df, ymin=min_df$Path, ymax=max_df$Path)
  # Create dataframe for cutoff
  cutoff <- data.frame(values=true_values, `Dimension`=c("X1", "X2", "X3"), stringsAsFactors = FALSE)
  # Plot
  ggplot(data=final_df) + 
    geom_ribbon(aes(x=index, ymin=ymin, ymax=ymax), fill="grey70") + 
    geom_line(aes(x=index, y=Path)) + 
    facet_wrap(~ `Dimension`, nrow=1) +
    geom_hline(data=cutoff, aes(yintercept=values), lty=2, alpha=0.5) +
    ggtitle("Path of Each Coordinate of the Posterior Mean") +
    theme(plot.title=element_text(hjust=0.5))
}


theta_posterior_path(true_values)