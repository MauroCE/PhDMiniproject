library(tidyverse)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(emdbook)
library(ellipse)
library(gridExtra)
library(ggpubr)


# UTILITY FUNCTIONS

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

natural <- function(mu, Sigma){    # Transforms (mu, Sigma) into a single natural parameter 
  Q <- matrix_inverse(Sigma)         # Could use chol2inv(chol()) but requires positive definiteness, and it breaks
  return(rbind(Q %*% mu, matrix(-0.5*Q, ncol=1)))   # TODO: Technically t(Q) but Q should be symmetric anyways?
}

natural2mean <- function(theta, d, returnMuSigma=FALSE){
  Sigma <- matrix_inverse(-2*matrix(theta[(d+1):(d+d^2), ], nrow=d))
  mu <- Sigma %*% matrix(theta[1:d, ])
  if (!returnMuSigma){
    return(moment(mu, Sigma))
  } else {
    return(list(mu=mu, Sigma=Sigma)) 
  }
}
to_exp <- function(theta, d){
  # transform to mu and Sigma
  muSigma <- natural2mean(theta, d, returnMuSigma = TRUE)
  # Change to log
  mu_exp    <- c(exp(muSigma$mu + 0.5*diag(muSigma$Sigma)))
  sigma_exp <- diag(mu_exp) %*% (exp(muSigma$Sigma) - 1) %*% diag(mu_exp)
  return(list(mu=matrix(mu_exp), Sigma=sigma_exp))
}
matrix_inverse <- function(matrix){
  m <- chol2inv(chol(matrix))
  if (class(m) != "matrix") return(pseudoinverse(matrix))
  else m
}
# FUNCTION TO READ IN A FILE
read_results <- function(file_name, path) unname(as.matrix(read.table(paste(path, file_name, sep=""))))
# FUNCTION TO TRANSFORM EACH LINE OF TP4W, ...
trans <- function(tt){
  n <- nrow(tt)
  means  <- matrix(0.0, nrow=n, ncol=3)
  sigmas    <- array(0.0, dim=c(3, 3, n))
  for (j in 1:n){
    ooo <- to_exp(matrix(tt[j, ]), 3)
    means[j, ] <- ooo$mu
    sigmas[,,j] <- ooo$Sigma
  }
  return(list(m=means, s=sigmas))
}

file2data <- function(filename, path){
  trans(read_results(filename, path))
}
# alternative
#toEllipse <- function(mu, Sigma, i) mutate(data.frame(ellipse::ellipse(cov2cor(Sigma), scale=sqrt(diag(Sigma)), center=mu)), iter=i)
toEllipse <- function(mu, Sigma, i) mutate(data.frame(ellipse(Sigma, centre=mu, level=0.95)), Iteration=i)

plotEvolution <- function(w, nworkers=2, dims=c(1,2), mu_prior=matrix(0, nrow=3), Sigma_prior=diag(3),
                          xrange=c(-5, 40), yrange=c(-5, 40), color_start="dodgerblue1",
                          color_final="darkblue", legend=TRUE){
  n <- nrow(w$m)
  # Data frame containing means (m) and Sigmas (s) of the evolution
  w_ellipse <- plyr::ldply(1:n, function(i) toEllipse(w$m[i, dims], w$s[dims,dims,i], i))
  # Grab mu and Sigma from start
  muexp_start <- to_exp(natural(mu_prior, Sigma_prior) * 3, 3)$mu
  Sigmaexp_start <- to_exp(natural(mu_prior, Sigma_prior) * 3, 3)$Sigma
  # Dataframe containing final solution
  #final <- toEllipse(w$m[n, dims], w$s[dims, dims, n], n)
  colored_df <- rbind(toEllipse(muexp_start, Sigmaexp_start, "Initial MVN"), 
                      toEllipse(w$m[n, dims], w$s[dims, dims, n], "Final MVN"))
  points <- data.frame(x=c(20, w$m[n, dims][1]), y=c(20, w$m[n, dims][2]), Points=c("Solution", "Final Mean"))
  p <- ggplot(data=w_ellipse) + 
    geom_path(aes(x, y, group=Iteration, alpha=Iteration), lwd=0.2) + 
    #geom_path(data=toEllipse(muexp_start, Sigmaexp_start, 0), aes(x, y), color=color_start) + 
    geom_path(data=colored_df, aes(x, y, color=Iteration), lwd=1) + #lwd=1, color=color_final) + 
    geom_point(data=points, aes(x=x, y=y, shape=Points)) +
    theme_minimal() +
    scale_x_continuous(limits=xrange) + 
    scale_y_continuous(limits=yrange) + 
    labs(x=expression(theta), y=expression(r)) +
    ggtitle(paste("Evolution with ", nworkers, " workers", sep="")) + 
    theme(plot.title=element_text(hjust=0.5)) + 
    coord_equal()
  if (!legend) p + theme(legend.position="none")
  else p #+ theme(legend.position="bottom")
}

#tp10w <- read_results("theta_posterior_history", path10w)
w2 <- file2data("theta_posterior_history", "~/Documents/University/PhDMiniproject/sfc_results_1_may/sfc_g_2w_40s_e005/")
w4 <- file2data("theta_posterior_history", "~/Documents/University/PhDMiniproject/sfc_results_1_may/sfc_g_4w_40s_e005/")
w8 <- file2data("theta_posterior_history", "~/Documents/University/PhDMiniproject/sfc_results_1_may/sfc_g_8w_40s_e005/")

p2 <- plotEvolution(w2, 2)
p2_legend <- get_legend(p2)
p2_new <- plotEvolution(w2, 2, legend=F)
p4 <- plotEvolution(w4, 4, legend=F)
p8 <- plotEvolution(w8, 8, legend=T)

#grid.arrange(arrangeGrob(p2_new, p4, p8, nrow=1), p2_legend, ncol=1) # + theme(legend.position = "bottom")
#grid.arrange(arrangeGrob(p2_new, p4, p8, nrow=1))
#grid.arrange(p2_new, p4, p8, nrow=1)

ggarrange(plotEvolution(w2, 2), 
          plotEvolution(w4, 4), 
          plotEvolution(w8, 8), 
          plotEvolution(w8, 8), 
          ncol=2, nrow=2, common.legend = TRUE, legend = "bottom") 





# 
# # I TAKE THE FIRST TWO DIMENSIONS
# dims <- 1:2
# #w2_ellipse <- ldply(1:nrow(w2$m), function(i) toEllipse(w2$m[i, ], w2$s[,,i], i))
# w2_ellipse <- plyr::ldply(1:33, function(i) toEllipse(w2$m[i, dims], w2$s[dims,dims,i], i))
# 
# mu_prior <- matrix(0, nrow=3)
# Sigma_prior <- diag(3)
# exp_stuff <- to_exp(natural(mu_prior, Sigma_prior) * 3, 3)
# mu_exp_start <- exp_stuff$mu
# Sigma_exp_start <- exp_stuff$Sigma
# final <- toEllipse(w2$m[33, dims], w2$s[dims, dims, 33], 33)
# points <- data.frame(x=c(20, w2$m[33, dims][1]), y=c(20, w2$m[33, dims][2]), Points=c("Solution", "Final Mean"))
# 
# ggplot(data=w2_ellipse) + 
#   geom_path(aes(x, y, group=Iteration, alpha=Iteration), lwd=0.2) + 
#   geom_path(data=toEllipse(mu_exp_start, Sigma_exp_start, 0), aes(x, y), color="dodgerblue1") + 
#   geom_path(data=final, aes(x, y), lwd=1, color="darkblue") + 
#   geom_point(data=points, aes(x=x, y=y, shape=Points), color=c("black", "darkblue")) +
#   theme_minimal() +
#   scale_x_continuous(limits=c(-5, 40)) + 
#   scale_y_continuous(limits=c(-5, 40)) + 
#   labs(x=expression(theta), y=expression(r))
# 
# 
# 
# # APPLY FUNCTION TO EACH THETA POSTERIOR
# ms2w <- trans(tp2w)
# ms4w <- trans(tp4w)
# ms8w <- trans(tp8w)
# #ms10w <- trans(tp10w)
# 
# # GRAB MEANS AND MAKE THEM DATAFRAMES
# todf <- function(m, nw){
#   out <- data.frame(workers=as.factor(nw), theta=m$m[, 1], r=m$m[, 2], "t[f]"=m$m[, 3])
#   names(out) <- c("workers", "theta", "r", "t[f]")
#   return(out %>% mutate(index=row_number()) %>% gather("coordinate", "value", -workers, -index))
# }
# dfm2 <- todf(ms2w, 2)
# dfm4 <- todf(ms4w, 4)
# dfm8 <- todf(ms8w, 8)
# #dfm10 <- todf(ms10w, 10)
# 
# # STACK DATAFRAMES TOGETHER
# df <- mutate(rbind(dfm4, dfm6, dfm8, dfm10))
# names(df) <- c("Workers", "Server Synchronizations", "Coordinate", "Value")
# 
# # PLOT
# ggplot(data=df, aes(x=`Server Synchronizations`, y=Value)) + 
#   geom_line(aes(linetype=Workers)) +
#   facet_wrap(~Coordinate, scales="free", labeller = label_parsed) + 
#   theme_minimal() + 
#   labs(title="Effect of Number of Workers on Convergence") + 
#   theme(strip.text.x = element_text(size=12),
#         plot.title = element_text(hjust=0.5))
# 
# 
# 
# # PLOT A CONTOUR PLOT FOR THE END RESULT
# m <- ms4w$m[nrow(ms4w$m), ]
# s <- diag(ms4w$s[nrow(ms4w$m), ]^2)
# # keep 1 and2 for now
# m <- m[c(1, 2)]
# s <- s[1:2, 1:2]
# # sample from it
# samples <- data.frame(mvrnorm(n=10000, mu=m, Sigma=s))
# data <- data.frame(samples, z=apply(samples, 1, function(x) dmvnorm(x, mu=m, Sigma=s)))
# ggplot(data=data, aes(x=X1, y=X2, z=z)) + geom_density_2d()
# 
# 
# 
# 
# # FIRST 4w
# m <- m
# sigma <- s
# data.grid <- expand.grid(s.1 = seq(0, 40, length.out=200), s.2 = seq(0, 40, length.out=200))
# q.samp <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m, sigma = sigma))
# # SECOND 6W
# m2 <- ms6w$m[nrow(ms6w$m), ][c(1, 2)]
# sigma2 <- diag(ms6w$s[nrow(ms6w$m), ]^2)[1:2, 1:2]
# q.samp2 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m2, sigma = sigma2))
# # third 8w
# m3 <- ms8w$m[nrow(ms8w$m), ][c(1, 2)]
# sigma3 <- diag(ms8w$s[nrow(ms8w$m), ]^2)[1:2, 1:2]
# q.samp3 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m3, sigma = sigma3))
# # fourth 10w
# m4 <- ms10w$m[nrow(ms10w$m), ][c(1, 2)]
# sigma4 <- diag(ms10w$s[nrow(ms10w$m), ]^2)[1:2, 1:2]
# q.samp4 <- cbind(data.grid, prob = mvtnorm::dmvnorm(data.grid, mean = m4, sigma = sigma4))
# 
# ggplot(q.samp, aes(x=s.1, y=s.2, z=prob)) + 
#   geom_contour(lty=1, bins=3) + 
#   geom_contour(data=q.samp2, aes(x=s.1, y=s.2, z=prob), lty=2, bins=3) +
#   geom_contour(data=q.samp3, aes(x=s.1, y=s.2, z=prob), lty=3, bins=3) + 
#   geom_contour(data=q.samp4, aes(x=s.1, y=s.2, z=prob), lty=4, bins=3) + 
#   coord_fixed(xlim=c(0,40), ylim=c(0,40),ratio=1)
# 
# 
# 
# # ALL OF THIS TOGETHER
# m4 <- ms4w$m
# S4 <- array(0.0, dim=c(2, 2, nrow(m4)))
# for (row_index in 1:nrow(m4)){
#   S4[,,row_index] <- diag(ms4w$s[row_index, 1:2]^2)
# }
# df_list <- list()
# alpha_levels <- seq(0.001, 0.95, by=0.8)
# names(alpha_levels) <- alpha_levels
# for (row_index in 1:nrow(m4)){
#   df_list[[row_index]] <- mutate(ldply(alpha_levels,ellipse,x=S4[,,row_index],
#                                 scale=c(1,1),  ## needed for positional matching
#                                 centre=m4[row_index, ]), iteration=row_index)
# }
# 
# ggplot() + 
#   geom_path(data=df_list[[1]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[2]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[3]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[4]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[5]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[6]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[7]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[8]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[9]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[10]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[11]], aes(x,y,group=.id, color=iteration)) + 
#   geom_path(data=df_list[[12]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[13]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[14]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[15]], aes(x,y,group=.id, color=iteration)) +
#   geom_path(data=df_list[[16]], aes(x,y,group=.id, color=iteration)) + 
#   theme_minimal() + 
#   labs(x=expression(theta), y=expression(r))
#   
# 
# 
# 
# 
# df <- data.frame(ellipse(cov2cor(Sigma), scale=sqrt(diag(Sigma)), centre=c(1, 2), level=lvl))
# ggplot() + 
#   geom_path(data=df, aes(x, y))  + geom_point(data=data.frame(x=1, y=2), aes(x, y))
# 
# 
# 
# 
# 
