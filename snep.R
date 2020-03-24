# MOMENT TO NATURAL PARAMETERS
moment2natural <- function(mu, Sigma){
  Q <- solve(Sigma)  # chol2inv(chol(Sigma)) is faster but sometimes fails
  r <- Q %*% mu
  return(list(r=r, Q=Q))
}
# Compute sufficient statistics (no need of t(tcrossprod(x)) because xx^t is symmetric)
suff <- function(x) rbind(x, matrix(tcrossprod(x), ncol=1))
# Find natural parameters from r and Q. No need to transpose Q as it should be symmetric anyways, but will keep it for now
natural <- function(r, Q) rbind(r, matrix(-0.5*t(Q), ncol=1))  # TODO: Does -0.5 matter here?
# Find mean parameters form mu and Sigma. No need to transpose second term as should be symmetric
moment <- function(mu, Sigma) rbind(mu, matrix(t(Sigma + tcrossprod(mu)), ncol=1))

natural2moment <- function(r, Q){
  Sigma <- solve(Q) # chol2inv(chol(Sigma)) is faster but sometimes fails
  mu <- Sigma %*% r
  return(list(mu=mu, Sigma=Sigma))
}

natural2moment_complete <- function(natural, d){
  # Separate natural into r and Q. Use dimensionality d
  r <- matrix(natural[1:d, ])
  Q <- natural[(d+1):(d+d^2), ]
  # Transform r and Q into moment parameter
  moment_params <- natural2moment(r, Q)
  # Transform mu and Sigma to full moment parameter
  return(moment(moment_params$mu, moment_params$Sigma))
}

moment2natural_complete <- function(moment_param, d){
  # Separate it into mu and Sigma
  mu <- matrix(moment_param[1:d, ])
  Sigma <- moment_param[(d+1):(d+d^2), ]
  # Find natural parameters from mu and Sigma, i.e. r and Q
  natural_param <- moment2natural(mu, Sigma)
  # Combine natural parameters
  return(natural(natural_param$r, natural_param$Q))
}

# RANDOM WAK METROPOLIS HASTINGS
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

######################################################################################################################################
######################################################################################################################################
####################################  A C T U A L    A L G O R I T H M ###############################################################
######################################################################################################################################
######################################################################################################################################
# SNEP MAIN ALGORITHM
snep <- function(nnodes, nouter, nsync, mu_local, Sigma_local, mu_prior, Sigma_prior,
                 sampler, betas, local_likelihood){
  # There must be a beta for every node
  if (length(betas) != n) stop("The number of beta_i's is not the same as the number of nodes.")
  
  # I N I T I A L I Z A T I O N
  # Convert local moment parameters to (initial) natural parameters
  natural_local <- moment2natural(mu_local, Sigma_local)
  r_local <- natural_local$r
  Q_local <- natural_local$Q
  # Convert prior moment parameters to (initial) natural parameters
  natural_prior <- moment2natural(mu_prior, Sigma_prior)
  r_prior <- natural_prior$r
  Q_prior <- natural_prior$Q
  # Find Initial natural parameters of cavity distribution. Use nnodes - 1 because at the beginning they are all the same
  r_c <- r_prior + (nnodes - 1) * r_local
  Q_c <- Q_prior + (nnodes - 1) * Q_local
  # Auxiliary Parameters should be the same at the beginning
  r_aux <- r_c + r_local
  Q_aux <- Q_c + Q_local
  # Sample a starting point from the proposal distribution with parameter equal to thetaprime (i.e. r_aux and Q_aux)
  
  # S E N D  S T U F F  T O  W O R K E R S
}


snep <- function(nnodes, nouter, nsync, mu_prior, Sigma_prior, mu_local, Sigma_local, betas, local_likelihood, t){
  # Grab dimensionality of parameter
  d <- nrow(mu_prior)
  # Prior: Natural Parameters theta_0
  rq_prior <- moment2natural(mu_prior, Sigma_prior)
  natural_prior <- natural(rq_prior$r, rq_prior$Q)    # theta_0
  # Local Approximation: Natural Parameters  
  rq_local <- moment2natural(mu_local, Sigma_local)
  natural_local <- natural(rq_local$r, rq_local$Q)    # lambda_i^(1)
  # Local Approximation: Mean Parameters
  moment_local <- moment(mu_local, Sigma_local)       # gamma_i^(1)
  # Since all sites are the same on initialization, we calculate cavity distribution parameters here? use initialize_workers?
  natural_cavity <- natural_prior + (nnodes - 1)*natural_local   # theta_{-i}
  # Find Auxiliary Parameter 
  aux <- natural_cavity + natural_local               # theta_i'
  # Sample initial starting point from a multivariate normal distribution with mean and variance found from the auxiliary 
  # natural parameter
  Sigma_aux <- -2*chol2inv(chol(matrix(aux[(d+1):(d^2+d)], nrow=d)))
  mu_aux    <- Sigma_aux %*% matrix(aux[1:d, ], ncol=1)
  initial_samples <- mvrnorm(n=nnodes, mu = mu_aux, Sigma = Sigma_aux)
  # Initialize workers
  workers <- list()
  for (i in nnodes){
    workers[[i]] <- initialize_worker(
      x_i=matrix(initial_samples[i, ]),
      gamma_i=moment_local,
      lambda_i=natural_local,
      aux_i = aux,
      beta_i = betas[i],
      natural_cavity=natural_cavity,
      t=t)
  }
  # Run workers asyncronously
  
} 


initialize_worker <- function(x_i, gamma_i, lambda_i, aux_i,
                              beta_i, natural_cavity, t){
  # Initialize worker
  worker <- list(
    x_i            = x_i,
    gamma_i        = gamma_i,
    lambda_i       = lambda_i,
    aux_i          = aux_i,
    beta_i         = beta_i,
    natural_cavity = natural_cavity,
    iteration      = 1,                    # stores iteration number
    t              = t                     # number of total iterations (rather than convergence)
    )
  return(worker)
}

# Use either i or worker directly
mainloop_worker <- function(i, t){
  # i : Index of the worker
  # t : Number of iterations for the main loop
  # Grab current worker
  w <- workers[[i]]
  # Sample x_i^(t+1) with RWMH targeting the log tilted
  log_tilted <- function(x)  (w$aux_i - (w$lambda_i / w$beta_i)) %*% suff(x) + (log(local_likelihood(x)) / x$beta_i)
  x_i_new <- rwmh(start = w$x_i, niter = 1, logtarget = log_tilted)
  # Use new sample to update gamma_i, the mean parameter of local approximation
  # TODO: either use 1/iteration or use something like 0.02, this is what they used in the paper
  w$gamma_i <- w$gamma_i + (1/w$iteration)*(suff(x_i_new) - natural2moment_complete(w$aux, nrow(w$aux)))
  # Compute natural parameters from the new gamma_i and update the worker
  w$lambda_i <- moment2natural_complete(w$gamma_i)
}




# worker <- function(mu_local, Sigma_local, r_local, Q_local, r_c, Q_c, beta_i, r_aux, Q_aux, xstart, local_likelihood){
#   # Log likelihood
#   log_likelihood <- function(x) log(local_likelihood(x))
#   # Construct Local Tilted Distribution to sample from
#   log_tilted <- function(x){
#     # calculate first parameter
#     r_first <- r_aux - r_local/beta_i
#     Q_first <- Q_aux - Q_local/beta_i
#     # Join parameters together
#     natural_first <- natural(r_first, Q_first)  # theta' - lambda_i^t / beta_i
#     # compute value
#     return(t(natural_first) %*% suff(x) + log_likelihood(x) / beta_i)
#   }
#   # Sample from tilted using RWMH
#   sample <- matrix(rwmh(start=xstart, niter=1, logtarget=log_tilted), ncol=1)
#   # compute sufficient statistics of this sample
#   
# }
