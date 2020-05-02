# FULL-COVARIANCE GAUSSIAN CASE
library(MASS)
#library(expm)
#library(pracma)
library(corpcor)

snep <- function(nworkers, nouter, nsync, mu_prior, Sigma_prior, mu_local, Sigma_local, betas, loglike,
                 epsilon, pop="stable", tol=1e-6, maxiter=1000, nsamples=1, digits=2, nfiles=100)
{
  # Dimensionality of parameter 
  d <- nrow(Sigma_prior)
  # Steps: 2 to 5 Initialize Natural and Mean parameters
  theta_0  <- natural(mu_prior, Sigma_prior)   
  lambda   <- natural(mu_local, Sigma_local)    
  gamma  <- moment(mu_local, Sigma_local)    
  theta_c  <- theta_0 + (nworkers - 1)*lambda      
  theta_p  <- theta_c + lambda                   
  # Sample start from a MVN with mean, vcov given by auxiliary parameter
  Sigma_p <- matrix_inverse(-2*matrix(theta_p[(d+1):(d^2+d), ], nrow=d))
  mu_p    <- Sigma_p %*% matrix(theta_p[1:d, ], ncol=1)
  starts <- mvrnorm(n=nworkers, mu = mu_p, Sigma = Sigma_p)
  # Initialize Theta posterior
  theta_posterior <-  theta_c + lambda
  # Store history of theta_posterior mean
  theta_posterior_history      <- matrix(0.0, nrow=((maxiter %/% nsync)+1), ncol=(d+d^2))
  theta_posterior_history[1, ] <- theta_posterior
  # Find site workers (chunks need to be calculated)
  chunken_sites <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE)) 
  chunk_sites   <- chunken_sites(1:nfiles, nworkers)
  # Store theta posterior for polyak averaging
  tp_polyak              <- theta_posterior
  tp_polyak_history      <- matrix(0.0, nrow=((maxiter %/% nsync)+1), ncol=(d+d^2))
  tp_polyak_history[1, ] <- tp_polyak
  # Store theta polyak average (that is, we keep a moving average of the theta posterior!)
  theta_posterior_polyak <- theta_posterior
  theta_posterior_polyak_history      <- matrix(0.0, nrow=((maxiter %/% nsync)+1), ncol=(d+d^2))
  theta_posterior_polyak_history[1, ] <- theta_posterior_polyak
  # Initialize workers
  workers <- list()
  for (i in 1:nworkers){
    # Store the starting x_i in the history of samples
    x_history      <- matrix(0.0, nrow=(maxiter+1), ncol=d)
    x_history[1, ] <- starts[i, ]
    # Store also the complete history. This contains all the samples determined by nsamples
    x_history_full      <- matrix(0.0, nrow=(maxiter*nsamples + 1), ncol=d)
    x_history_full[1, ] <- starts[i, ]
    # Store history of the local likelihood approximation lambda_i, averaged
    workers[[i]]   <- list(x              = matrix(starts[i, ]),
                           gamma          = gamma,
                           gamma_avg      = gamma,
                           lambda         = lambda,
                           lambda_old     = lambda,
                           theta_p        = theta_p,
                           beta           = betas[i],
                           theta_c        = theta_c,
                           delta          = rep(Inf, d+d^2),              # Initialize the delta at Inf
                           sites          = chunk_sites[[i]],
                           x_history      = x_history,
                           x_history_full = x_history_full,
                           write_h        = F,
                           write_h_full   = F
    )                
  }
  # Run workers SYNCHRONOUSLY
  iter     <- 1
  distance <- Inf
  # Debug
  cat("Iter: ", iter, " Sync: ", iter %/% nsync, 
      " Dist: ", format(distance, digits=digits), 
      " Mu Polyak: ", format(to_exp(tp_polyak, d)$mu, digits=digits),
      " Mu Post: ", format(to_exp(theta_posterior, d)$mu, digits=digits), 
      " Mu Post Polyak: ", format(to_exp(theta_posterior_polyak, d)$mu, digits=digits), "\n")
  while (distance > tol){
    # INNER LOOP (WORKERS LEARN)
    for (j in 1:nworkers){
      muSigma_c <- natural2mean(workers[[j]]$theta_c, d, returnMuSigma = TRUE)  # mu, Sigma of cavity distribution
      if (!is_positive_definite(muSigma_c$Sigma)){                   # Sigma_cavity must be positive definite
        cat("Skipping worker", j, "at iteration", iter, "as Sigma of cavity distribution is not positive definite.\n")
        cat(muSigma_c$Sigma, "\n")
        next
      }
      # Sample from log tilted distribution
      tilted_param <- workers[[j]]$theta_p - (workers[[j]]$lambda / workers[[j]]$beta)
      logtilted <- function(x){
        return(as.double(t(tilted_param) %*% s(x)) + loglike(j, workers[[j]]$sites, c(x), pop=pop) / workers[[j]]$beta)
      }   
      # Average multiple samples (drop=F required if nsamples=1)
      samples <- rwmh(start=workers[[j]]$x, niter=(nsamples+1), logtarget=logtilted)[2:(nsamples+1), ,drop=F]
      # store all samples into the full history
      workers[[j]]$x_history_full[(2+nsamples*(iter-1)):(1+nsamples*iter), ] <- samples
      # WRITE FULL HISTORY
      write(t(samples), file=paste("x_h_full_w", j, sep=""), ncolumns=d, append=workers[[j]]$write_h_full)
      workers[[j]]$write_h_full <- T
      # compute monte carlo average of the sufficient statistics
      avg_suff_x <- matrix(0.0, nrow=(d+d^2))
      for (ii in 1:nsamples){
        avg_suff_x <- avg_suff_x + s(matrix(samples[ii, ]))
      }
      avg_suff_x <- avg_suff_x / nsamples
      # Update x for worker and history of x
      workers[[j]]$x                   <- matrix(samples[nsamples, ])
      workers[[j]]$x_history[iter+1, ] <- samples[nsamples, ]
      write(samples[nsamples, ], file=paste("x_h_w", j, sep=""), ncolumns=d, append=workers[[j]]$write_h)
      workers[[j]]$write_h <- T
      # Update gamma and lambda (mean and natural param of local approximation respectively)
      workers[[j]]$gamma  <- workers[[j]]$gamma + epsilon*(avg_suff_x - natural2mean(workers[[j]]$theta_c + workers[[j]]$lambda, d))
      # Polyak averaging
      workers[[j]]$gamma_avg <- workers[[j]]$gamma_avg + (workers[[j]]$gamma - workers[[j]]$gamma_avg) / (iter + 1)
      workers[[j]]$lambda <- mean2natural(workers[[j]]$gamma, d)
    }
    # OUTER LOOP (UPDATE AUXILIARY)
    if (iter %% nouter == 0) {
      for (j in 1:nworkers){
        workers[[j]]$theta_p <- workers[[j]]$theta_c + workers[[j]]$lambda
      }
    }
    # SYNCHRONIZATION WITH SERVER (UPDATE CAVITY and THETA POSTERIOR)
    if (iter %% nsync  == 0){
      for (j in 1:nworkers){
        # Worker updates delta and lambda old. Sends delta to posterior
        workers[[j]]$delta          <- workers[[j]]$lambda - workers[[j]]$lambda_old  
        workers[[j]]$lambda_old     <- workers[[j]]$lambda                 
        # Server updates posterior
        theta_posterior             <- theta_posterior + workers[[j]]$delta
        # Worker receives back theta posterior and updates cavity distribution (prior)
        # Also shift the state of the MCMC to avoid using burn-in
        #old <- natural2mean(workers[[j]]$theta_c + workers[[j]]$lambda, d, returnMuSigma = TRUE)
        workers[[j]]$theta_c <- theta_posterior - workers[[j]]$lambda_old
        # only update if both matrices allow for cholesky decomposition
        #new <- natural2mean(workers[[j]]$theta_c + workers[[j]]$lambda, d, returnMuSigma = TRUE)
        #workers[[j]]$x <- new$mu + sqrtm(new$Sigma) %*% matrix_inverse(sqrtm(old$Sigma), workers[[j]]$x - old$mu)
      }
      # Update history of theta posterior
      theta_posterior_history[(iter %/% nsync)+1, ] <- theta_posterior
      # Update theta polyak average
      theta_posterior_polyak <- theta_posterior_polyak + (theta_posterior - theta_posterior_polyak) / ((iter %/% nsync) + 1)
      # Update history of theta posterior polyak (simple average of theta posterior)
      theta_posterior_polyak_history[(iter %/% nsync)+1, ] <- theta_posterior_polyak
      # Find new theta posterior from polyak averaging (the theoretically correct one)
      tp_polyak_new <- theta_0
      for (j in 1:nworkers) tp_polyak_new <- tp_polyak_new + mean2natural(workers[[j]]$gamma_avg, d)
      # Update distance
      distance <- sqrt(sum((tp_polyak - tp_polyak_new)^2))
      # Update theta posterior from polyak and its history
      tp_polyak                               <- tp_polyak_new
      tp_polyak_history[(iter %/% nsync)+1, ] <- tp_polyak
      # WRITE EVERYTHING TO FILES
      append_flag <- (nsync / iter) != 1
      write(theta_posterior, file="theta_posterior_history", ncolumns=(d+d^2), append=append_flag)
      write(theta_posterior_polyak, file="theta_posterior_polyak", ncolumns=(d+d^2), append=append_flag)
      write(tp_polyak, file="tp_polyak_history", ncolumns=(d+d^2), append=append_flag)
      # Debug
      cat("Iter: ", iter, " Sync: ", iter %/% nsync, 
          " Dist: ", format(distance, digits=digits), 
          " MuPol: ", format(to_exp(tp_polyak, d)$mu, digits=digits),
          " MuPost: ", format(to_exp(theta_posterior, d)$mu, digits=digits), 
          " MuPostPol: ", format(to_exp(theta_posterior_polyak, d)$mu, digits=digits), "\n")
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
              theta_posterior_history=theta_posterior_history,
              tp_polyak_history=tp_polyak_history,
              theta_posterior_polyak=theta_posterior_polyak,
              theta_posterior_polyak_history=theta_posterior_polyak_history))
} 


### HELPER FUNCTIONS
matrix_inverse <- function(matrix){
  m <- chol2inv(chol(matrix))
  if (class(m) != "matrix") return(pseudoinverse(matrix))
  else m
}

natural <- function(mu, Sigma){    # Transforms (mu, Sigma) into a single natural parameter 
  Q <- matrix_inverse(Sigma)         # Could use chol2inv(chol()) but requires positive definiteness, and it breaks
  return(rbind(Q %*% mu, matrix(-0.5*Q, ncol=1)))   # TODO: Technically t(Q) but Q should be symmetric anyways?
}
moment <- function(mu, Sigma){     # Transforms (mu, Sigma) into the mean parameters for a Gaussian distribution
  return(rbind(mu, matrix(Sigma + tcrossprod(mu), ncol=1))) # TODO: Technically t(Sigma + tcrossprod(mu)) but symmetric
}
natural2mean <- function(theta, d, returnMuSigma=FALSE){   # Transforms a complete natural parameter to mean parameter for Gaussian Approx
  #Sigma <- -0.5*matrix_inverse(-2*matrix(theta[(d+1):(d+d^2), ], nrow=d))  # as (-2)^{-1} = -0.5
  Sigma <- matrix_inverse(-2*matrix(theta[(d+1):(d+d^2), ], nrow=d))
  mu <- Sigma %*% matrix(theta[1:d, ])
  if (!returnMuSigma){
    return(moment(mu, Sigma))
  } else {
    return(list(mu=mu, Sigma=Sigma)) 
  }
}
mean2natural <- function(nu, d){  # Given a (d + d^2) dimensional mean param, it returns its correponding natural param
  mu    <- matrix(nu[1:d, ])
  Sigma <- matrix(nu[(d+1):(d+d^2), ], nrow=d) - tcrossprod(mu)
  return(natural(mu, Sigma))
}
s <- function(x) rbind(matrix(x), matrix(tcrossprod(matrix(x)), ncol=1)) # Calculate sufficient statistics
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

to_exp <- function(theta, d){
  # transform to mu and Sigma
  muSigma <- natural2mean(theta, d, returnMuSigma = TRUE)
  # Change to log
  mu_exp    <- c(exp(muSigma$mu + 0.5*diag(muSigma$Sigma)))
  sigma_exp <- diag(mu_exp) %*% (exp(muSigma$Sigma) - 1) %*% diag(mu_exp)
  return(list(mu=matrix(mu_exp), Sigma=sigma_exp))
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

## LOG LIKELIHOOD
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
nworkers    <- 4                                  # Number of sites
nouter      <- 3                                 # Number of iterations after which we update auxiliary parameter
nsync       <- 3                                 # Number of iterations after which we communicate with the server
mu_prior    <- mu_local    <- matrix(rep(0, d))
Sigma_prior <- Sigma_local <- diag(d)
betas       <- rep(1, nworkers)          # pSNEP
epsilon     <- 0.05                               # learning rate
pop         <- "growing"                      # Choose between "stable", "growing", or "contracting"
tol         <- 1e-6                               # Tolerance used to determined whether theta_posterior has converged
maxiter     <- 100
nsamples    <- 40
digits      <- 4
nfiles      <- 100
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
            maxiter=maxiter,
            nsamples=nsamples,
            digits=digits,
            nfiles=nfiles
)



