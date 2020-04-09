# DIAGONAL GAUSSIAN CASE
library(MASS)
# USING HIS SINTAX
snep <- function(nworkers, nouter, nsync, mu_prior, Sigma_prior, mu_local, Sigma_local, betas, loglike,
                 epsilon, pop="stable", tol=1e-6, maxiter=1000)
{
  # Grab dimension of parameter space (usually 3)
  d       <- nrow(Sigma_prior)
  # Initialize Natural and Mean Parameters
  theta_0 <- natural(mu_prior, Sigma_prior)
  lambda  <- natural(mu_local, Sigma_local)
  gamma   <- moment(lambda, d)
  theta_c <- theta_0 + (nworkers-1)*lambda
  theta_p <- theta_c + lambda
  # Sample initial states of MCMC samplers
  starts <- mvrnorm(n = nworkers, mu = -theta_p[1:d]/theta_p[(d+1):(2*d)], Sigma = diag(-1/theta_p[(d+1):(2*d)]))
  # Initialize Theta Posterior
  theta_posterior <- theta_c + lambda
  # Story history of theta posterior
  theta_posterior_history      <- matrix(0.0, nrow=(maxiter %/% nsync+1), ncol=(2*d))
  theta_posterior_history[1, ] <- theta_posterior
  # Find site workers
  chunken_sites <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  chunk_sites   <- chunken_sites(1:20, nworkers)
  # store thera posterior for polyak averaging
  tp_polyak         <- theta_posterior
  tp_polyak_history      <- matrix(0.0, nrow=(maxiter %/% nsync+1), ncol=(2*d))
  tp_polyak_history[1, ] <- theta_posterior
  # Initialize Workers
  workers <- list()
  for (i in 1:nworkers){
    # Store the starting x_i in the history of samples
    x_history      <- matrix(0.0, nrow=(maxiter+1), ncol=d)
    x_history[1, ] <- starts[i, ]
    # Store history of the local likelihood approximation lambda_i, averaged. 
    workers[[i]]   <- list(x          = matrix(starts[i, ]),
                           gamma      = gamma,
                           gamma_avg  = gamma,
                           lambda     = lambda,
                           lambda_old = lambda,
                           theta_p    = theta_p,
                           beta       = betas[i],
                           theta_c    = theta_c,
                           delta      = rep(Inf, 2*d),              # Initialize the delta at Inf
                           sites      = chunk_sites[[i]],             # Indeces of site chunks
                           x_history  = x_history
    )                
  }
  # Run workers SYNCHRONOUSLY
  iter     <- 1
  distance <- Inf
  # Debug
  cat("Iteration: ", iter, " Theta_Update: ", iter %/% nsync, 
      " Distance: ", distance, " Mu: ", to_exp(tp_polyak, d)$mu, "\n")
  while (distance > tol){
    # INNER LOOP (WORKERS LEARN)
    for (j in 1:nworkers){
      # Sample from log tilted distribution
      tilted_param <- workers[[j]]$theta_p - workers[[j]]$lambda/workers[[j]]$beta
      log_tilted <- function(x){  # TODO: Check if I need to do c(x)
        return(as.double(tilted_param %*% s(x)) + loglike(j, workers[[j]]$sites, x, pop=pop) / workers[[j]]$beta)
      }
      x_new <- rwmh(start=workers[[j]]$x, niter=2, logtarget=log_tilted)[2, ]
      # Update x for worker and history of x
      workers[[j]]$x                   <- x_new
      workers[[j]]$x_history[iter+1, ] <- x_new
      # Update natural and mean parameters
      workers[[j]]$gamma  <- workers[[j]]$gamma + (iter^(-2/3))*(s(x_new) - moment(workers[[j]]$theta_c + workers[[j]]$lambda, d))
      # use polyak averaging
      workers[[j]]$gamma_avg <- workers[[j]]$gamma_avg + (workers[[j]]$gamma - workers[[j]]$gamma_avg) / (iter + 1)
      workers[[j]]$lambda <- moment2natural(workers[[j]]$gamma, d)
    }
    # OUTER LOOP  (UPDATE AUXILIARY)
    if (iter %% nouter == 0){
      for (j in 1:nworkers){
        workers[[j]]$theta_p <- workers[[j]]$theta_c + workers[[j]]$lambda
      }
    }
    # SYNCHRONIZATION WITH SERVER (UPDATE CAVITY and THETA POSTERIOR)
    if (iter %% nsync == 0){
      for (j in 1:nworkers){
        # Worker updates delta and lambda old. Sends delta to posterior
        workers[[j]]$delta      <- workers[[j]]$lambda - workers[[j]]$lambda_old
        workers[[j]]$lambda_old <- workers[[j]]$lambda
        # Server updates posterior
        theta_posterior         <- theta_posterior + workers[[j]]$delta
        # Worker receives back theta posterior and updates cavity distribution (prior)
        # Also shift the state of the MCMC to avoid using burn-in
        old <- natural2muSigma(workers[[j]]$theta_c + workers[[j]]$lambda, d)
        workers[[j]]$theta_c <- theta_posterior - workers[[j]]$lambda_old         # UPDATE CAVITY
        new <- natural2muSigma(workers[[j]]$theta_c + workers[[j]]$lambda, d)
        workers[[j]]$x <- new$mu + sqrt(new$sigmas)*sqrt(1/old$sigmas)*(workers[[j]]$x - old$mu)   # SHIFT MCMC STATE
      }
      # Update history of theta posterior
      theta_posterior_history[(iter %/% nsync)+1, ] <- theta_posterior
      # Find new theta posterior from polyak averaging
      tp_polyak_new <- theta_0
      for (j in 1:nworkers) tp_polyak_new <- tp_polyak_new + moment2natural(workers[[j]]$gamma_avg, d)
      # Update distance
      distance <- sqrt(sum((tp_polyak - tp_polyak_new)^2))
      # Update theta posterior polyak, and its history
      tp_polyak <- tp_polyak_new
      tp_polyak_history[(iter %/% nsync)+1, ] <- tp_polyak
      # Debug
      cat("Iteration: ", iter, " Theta_Update: ", iter %/% nsync, 
          " Distance: ", distance, " Mu: ", to_exp(tp_polyak, d)$mu, "\n")
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
  return(list(workers=workers, server=theta_posterior, 
              theta_posterior_history=theta_posterior_history, tp_polyak_history=tp_polyak_history))
}
# Function to take theta posterior and (or any natural parameter) and transform it into mean and variance covariance
# martrix of phi
to_exp <- function(theta, d){
  # transform to mu and Sigma
  muSigma <- natural2muSigma(theta, d)
  # Change to log
  mu_exp    <- exp(muSigma$mu + 0.5*muSigma$sigmas)
  sigma_exp <- diag(mu_exp) %*% (exp(diag(muSigma$sigmas)) - 1) %*% diag(mu_exp)
  return(list(mu=mu_exp, Sigma=sigma_exp))
}



# UTILITIES
# Natural Parameters from mu, Sigma
natural <- function(mu, Sigma){
  return(c(mu/diag(Sigma), -1/diag(Sigma)))
}
# Mean Parameters from Natural Parameters
moment <- function(theta, d){
  mu <- -theta[1:d]/theta[(d+1):(2*d)]
  return(c(mu, 0.5*(mu*mu - 1/theta[(d+1):(2*d)])))
}
# From moment parameter to natural parameter
moment2natural <- function(eta, d){
  scaling <- 1 / (2*eta[(d+1):(2*d)] - eta[1:d]^2)
  return(c(eta[1:d]*scaling, -scaling))
}
natural2muSigma <- function(theta, d){
  diagonal <- -1/theta[(d+1):(2*d)]
  return(list(mu=theta[1:d]*diagonal, sigmas=diagonal))
}
# From mu, Sigma to natural (used for sampling at the beginning)
s <- function(x) c(x, 0.5*x^2)

# SAMPLING
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

# LOG LIKELIHOOD
loglike <- function(worker_index, sites, params, isss=5000, pop="stable"){
  # Name of file containing phi (x)
  paramfile <- paste("paramfiles/paramfile_new", worker_index, sep="")
  write(params, paramfile, ncol=length(params)) 
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
nworkers    <- 2                                  # Number of sites
nouter      <- 10                                 # Number of iterations after which we update auxiliary parameter
nsync       <- 10                                 # Number of iterations after which we communicate with the server
mu_prior    <- mu_local    <- matrix(rep(0, d))
Sigma_prior <- Sigma_local <- diag(d)
betas       <- rep(1, nworkers, nworkers)          #rep(1, nworkers)                  #rep(1/nworkers, nworkers)          # pSNEP
epsilon     <- 0.05                               # learning rate
pop         <- "growing"                      # Choose between "stable", "growing", or "contracting"
tol         <- 1e-6                               # Tolerance used to determined whether theta_posterior has converged
maxiter     <- 1000                         
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


# write polyak posterior history
write(t(out$tp_polyak_history), file=paste("polyak_history_", pop, sep=""), ncolumns=6)
# write theta posterior history (not averaged)
write(t(out$theta_posterior_history), file=paste("tp_history_", pop, sep=""), ncolumns=6)
# write trace of every worker (here we have only 2 workers, so only write first two)
write(t(out$workers[[1]]$x_history), file=paste("x_history_w1_", pop, sep=""), ncolumns=3)
write(t(out$workers[[2]]$x_history), file=paste("x_history_w2_", pop, sep=""), ncolumns=3)











