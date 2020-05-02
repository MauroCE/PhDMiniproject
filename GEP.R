library(MASS)
library(pracma)

# Transforms (mu, Sigma) into (r, Q) where Q = Sigma^{-1} and r = Sigma^{-1} mu
natural <- function(mu, Sigma) {Q = matrix_inverse(Sigma); list(r =(Q %*% mu),    Q    =Q)}
muSigma <- function(r, Q)      {Sigma = matrix_inverse(Q); list(mu=(Sigma %*% r), Sigma=Sigma)}
suff    <- function(x) rbind(x, matrix(tcrossprod(x), ncol=1))


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

is_positive_definite <- function(Sigma, tol=1e-12){
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  bool <- !all(ev >= -tol*abs(ev[1L]))
  ifelse(bool, F, T)
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

### LOG LIKELIHOOD FUNCTION
loglike <- function(site, params, isss=5000, pop="stable"){
  # Name of file containing phi (x)
  paramfile <- paste("paramfiles_gep/paramfile_new", site, sep="")
  write(params, paramfile, ncol=length(params)) 
  # Create and store file names for target data points of chunk
  target <- paste(pop, "/gtree_file", site, sep="")
  # Now compute the product of all the likelihoods for the chunk
  ll <- log(as.numeric(system(paste("./is_moments_new",target,as.integer(1),as.integer(isss),paramfile, T), intern=T)))
  return(ll)
}

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
niter       <- 10
nsamples    <- 100
nfiles      <- 50
mu_prior    <- matrix(rep(0, 3))
Sigma_prior <- diag(3)
pop         <- "growing"
alpha       <- 0.01    # damping factor

################################################# RUN
### INITIALIZATION
# Transform to Q and r parameters
rQprior <- natural(mu_prior, Sigma_prior)
r_prior  <- rQprior$r 
Q_prior  <- rQprior$Q
# Store local parameters for each loci (100)
r_sites  <- array(0.0, dim=c(3, 1, nfiles))
Q_sites  <- array(0.0, dim=c(3, 3, nfiles))
# Form global approximation parameter (without the prior, it will be added at the end!)
r_global <- r_prior #apply(r_sites, c(1, 2), sum) + r_prior
Q_global <- Q_prior #apply(Q_sites, c(1, 2), sum) + Q_prior
# MAIN LOOP
for (i in 1:niter){
  cat("### Iteration: ", i, "\n")
  for (site in 1:nfiles){
    # Form cavity (r, Q)
    r_cavity <- r_global - matrix(r_sites[,,site])
    Q_cavity <- Q_global - Q_sites[,,site]
    # Transform to (mu, Sigma)
    muSigma_cavity <- muSigma(r_cavity, Q_cavity)
    mu_cavity      <- muSigma_cavity$mu
    Sigma_cavity   <- muSigma_cavity$Sigma
    # Sigma of cavity distribution needs to be positive definite. Otherwise, skip update
    if (!is_positive_definite(Sigma_cavity, tol = 1e-9)) {
      cat("Skipping site ", site, "\n")
      next
    }
    # Sample from cavity distribution to have a starting point for RWMH
    #start <- mvrnorm(n=1, mu=mu_cavity, Sigma=Sigma_cavity)
    # Form log-tilted distribution
    logtilted <- function(x){
      -0.5*t(x) %*% (Q_cavity %*% x) + t(r_cavity) %*% x + loglike(site, x, pop=pop)
    }
    # Sample from the tilted distribution
    samples <- rwmh(start=mu_cavity, niter=(nsamples+1), logtarget=logtilted, Sigma=diag(3))[2:(nsamples+1), ]
    # Average sufficient statistics
    #rQ_tilted <- matrix(0.0, nrow=)
    #for (jj in 1:nsamples) rQ_tilted <- rQ_tilted + suff(matrix(samples[jj, ], ncol=1))
    Q_tilted <- chol2inv(qr.R(qr(center_colmeans(samples)))) * (nsamples - 3 - 2)
    r_tilted <- Q_tilted %*% matrix(apply(samples, 2, mean))
    # Obtain empirical estimates of mu_tilted and Sigma_tilted
    #mu_tilted    <- apply(samples, 2, mean)
    #Sigma_tilted <- cov(samples)
    # Transform tilted (mu, Sigma) into (r, Q)
    #rQtilted <- natural(mu_tilted, Sigma_tilted)
    #r_tilted <- rQtilted$r
    #Q_tilted <- rQtilted$Q
    r_sites[,,site] <- (1-alpha)*r_sites[,,site] + alpha*(r_tilted - r_cavity)
    Q_sites[,,site] <- (1-alpha)*Q_sites[,,site] + alpha*(Q_tilted - Q_cavity)
    r_global <- apply(r_sites, c(1,2), sum) + r_prior
    Q_global <- apply(Q_sites, c(1,2), sum) + Q_prior
    # Update (r, Q) global
    #r_global <- alpha*r_tilted + (1-alpha)*r_global
    #Q_global <- alpha*Q_tilted + (1-alpha)*Q_global
    # Update site (r, Q)
    #r_sites[,,site] <- r_global - r_cavity
    #Q_sites[,,site] <- Q_global - Q_cavity
    # find parameters exp
    muSigmaExp <- to_exp2(r_global, Q_global)
    cat("Iter: ", i, "Site: ", site, "muExp: ", muSigmaExp$mu, "Sigma: ", diag(muSigmaExp$Sigma),"\n")
  }
  muSigmaExp <- to_exp2(r_global, Q_global)
  cat("mu: ", muSigmaExp$mu, "Sigma: ", diag(muSigmaExp$Sigma), "\n")
}
# Transform global (r, Q) to (mu, Sigma)
muSigma_global  <- muSigma(r_global, Q_global)
mu_global       <- muSigma_global$mu
Sigma_global    <- muSigma_global$Sigma


