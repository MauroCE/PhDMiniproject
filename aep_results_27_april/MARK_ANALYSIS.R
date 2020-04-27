#### ANALYZING MARK RESULTS

# Read Data
file_path <- "~/Documents/University/PhDMiniproject/aep_results_27_april/mark_contracting/contracting_mark.txt"
m <- unname(as.matrix(read.table(file_path)))

# For each iteration, build up muExp and SigmaExp
toExp <- function(i){
  # Grab mu and Sigma
  mu    <- matrix(m[i, 4:6])
  Sigma <- matrix(m[i, 7:15], nrow=3, ncol=3)
  # Find muExp and SigmaExp
  muExp    <- c(exp(mu + 0.5*diag(Sigma)))
  SigmaExp <- diag(muExp) %*% (exp(Sigma) - 1) %*% diag(muExp)
  # Combine them for ease
  combined <- c(muExp, c(SigmaExp))
  return(combined)
}


# Now apply function to m and extract mean and sigmas
applied <- lapply(1:nrow(m), toExp)
mus <- matrix(0.0, nrow=nrow(m), ncol=3)
sigmas <- array(0.0, dim = c(3, 3, nrow(m)))
for (iter in 1:nrow(m)){
  mus[iter, ] <- applied[[iter]][1:3]
  sigmas[,,iter] <- matrix(applied[[iter]][4:12], nrow=3, ncol=3)
}