library(tidyverse)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(emdbook)

# SETTINGS
path2w = paste("~/Documents/University/PhDMiniproject/snep_results_multiple_workers/sfc_g_2w_40s_e005//")
path4w = paste("~/Documents/University/PhDMiniproject/snep_results_multiple_workers/sfc_g_4w_40s_e005//")
path8w = paste("~/Documents/University/PhDMiniproject/snep_results_multiple_workers/sfc_g_8w_40s_e005//")
#path10w = paste("~/Documents/University/PhDMiniproject/snep_results_multiple_workers/growing_10w/")

# GET TRUE VALUES
true_values <- c(20, 20, 0.05) # stable c(10, 1, 1)   # contracting c(1, 0.05, 1)

# FUNCTION TO READ IN A FILE
read_results <- function(file_name, path) unname(as.matrix(read.table(paste(path, file_name, sep=""))))

# READ THETA POSTERIOR ONLY
tp2w <- read_results("theta_posterior_history", path2w)
tp4w <- read_results("theta_posterior_history", path4w)
tp8w <- read_results("theta_posterior_history", path8w)
#tp10w <- read_results("theta_posterior_history", path10w)

# FUNCTION TO TRANSFORM EACH LINE OF TP4W, ...
trans <- function(tt){
  n <- nrow(tt)
  means  <- matrix(0.0, nrow=n, ncol=3)
  sigmas    <- array(0.0, dim=c(3, 3, n))
  for (j in 1:n){
    mean_params <- natural2muSigma(matrix(tt[j, ]), 3)
    means[j, ] <- exp(mean_params$mu + 0.5*mean_params$sigmas)
    sigmas[,,j] <- diag(means[j, ]) %*% (exp(diag(mean_params$sigmas)) - 1) %*% diag(means[j, ])
  }
  return(list(m=means, s=sigmas))
}

# APPLY FUNCTION TO EACH THETA POSTERIOR
ms2w <- trans(tp2w)
ms4w <- trans(tp4w)
ms8w <- trans(tp8w)
#ms10w <- trans(tp10w)

# GRAB MEANS AND MAKE THEM DATAFRAMES
todf <- function(m, nw){
  out <- data.frame(workers=as.factor(nw), theta=m$m[, 1], r=m$m[, 2], "t[f]"=m$m[, 3])
  names(out) <- c("workers", "theta", "r", "t[f]")
  return(out %>% mutate(index=row_number()) %>% gather("coordinate", "value", -workers, -index))
}
dfm4 <- todf(ms4w, 4)
dfm6 <- todf(ms6w, 6)
dfm8 <- todf(ms8w, 8)
#dfm10 <- todf(ms10w, 10)

# STACK DATAFRAMES TOGETHER
df <- mutate(rbind(dfm4, dfm6, dfm8, dfm10))
names(df) <- c("Workers", "Server Synchronizations", "Coordinate", "Value")

# PLOT
ggplot(data=df, aes(x=`Server Synchronizations`, y=Value)) + 
  geom_line(aes(linetype=Workers)) +
  facet_wrap(~Coordinate, scales="free", labeller = label_parsed) + 
  theme_minimal() + 
  labs(title="Effect of Number of Workers on Convergence") + 
  theme(strip.text.x = element_text(size=12),
        plot.title = element_text(hjust=0.5))



# PLOT A CONTOUR PLOT FOR THE END RESULT
m <- ms4w$m[nrow(ms4w$m), ]
s <- diag(ms4w$s[nrow(ms4w$m), ]^2)
# keep 1 and2 for now
m <- m[c(1, 2)]
s <- s[1:2, 1:2]
# sample from it
samples <- data.frame(mvrnorm(n=10000, mu=m, Sigma=s))
data <- data.frame(samples, z=apply(samples, 1, function(x) dmvnorm(x, mu=m, Sigma=s)))
ggplot(data=data, aes(x=X1, y=X2, z=z)) + geom_density_2d()




# FIRST 4w
m <- m
sigma <- s
data.grid <- expand.grid(s.1 = seq(0, 40, length.out=200), s.2 = seq(0, 40, length.out=200))
q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))
# SECOND 6W
m2 <- ms6w$m[nrow(ms6w$m), ][c(1, 2)]
sigma2 <- diag(ms6w$s[nrow(ms6w$m), ]^2)[1:2, 1:2]
q.samp2 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m2, sigma = sigma2))
# third 8w
m3 <- ms8w$m[nrow(ms8w$m), ][c(1, 2)]
sigma3 <- diag(ms8w$s[nrow(ms8w$m), ]^2)[1:2, 1:2]
q.samp3 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m3, sigma = sigma3))
# fourth 10w
m4 <- ms10w$m[nrow(ms10w$m), ][c(1, 2)]
sigma4 <- diag(ms10w$s[nrow(ms10w$m), ]^2)[1:2, 1:2]
q.samp4 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m4, sigma = sigma4))

ggplot(q.samp, aes(x=s.1, y=s.2, z=prob)) + 
  geom_contour(lty=1, bins=3) + 
  geom_contour(data=q.samp2, aes(x=s.1, y=s.2, z=prob), lty=2, bins=3) +
  geom_contour(data=q.samp3, aes(x=s.1, y=s.2, z=prob), lty=3, bins=3) + 
  geom_contour(data=q.samp4, aes(x=s.1, y=s.2, z=prob), lty=4, bins=3) + 
  coord_fixed(xlim=c(0,40), ylim=c(0,40),ratio=1)



# ALL OF THIS TOGETHER
m4 <- ms4w$m
S4 <- array(0.0, dim=c(2, 2, nrow(m4)))
for (row_index in 1:nrow(m4)){
  S4[,,row_index] <- diag(ms4w$s[row_index, 1:2]^2)
}
df_list <- list()
alpha_levels <- seq(0.001, 0.95, by=0.8)
names(alpha_levels) <- alpha_levels
for (row_index in 1:nrow(m4)){
  df_list[[row_index]] <- mutate(ldply(alpha_levels,ellipse,x=S4[,,row_index],
                                scale=c(1,1),  ## needed for positional matching
                                centre=m4[row_index, ]), iteration=row_index)
}

ggplot() + 
  geom_path(data=df_list[[1]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[2]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[3]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[4]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[5]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[6]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[7]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[8]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[9]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[10]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[11]], aes(x,y,group=.id, color=iteration)) + 
  geom_path(data=df_list[[12]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[13]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[14]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[15]], aes(x,y,group=.id, color=iteration)) +
  geom_path(data=df_list[[16]], aes(x,y,group=.id, color=iteration)) + 
  theme_minimal() + 
  labs(x=expression(theta), y=expression(r))
  

# alternative
toEllipse <- function(mu, Sigma) data.frame(ellipse(cov2cor(Sigma), scale=sqrt(diag(Sigma)), center=mu))



df <- data.frame(ellipse(cov2cor(Sigma), scale=sqrt(diag(Sigma)), centre=c(1, 2), level=lvl))
ggplot() + 
  geom_path(data=df, aes(x, y))  + geom_point(data=data.frame(x=1, y=2), aes(x, y))





