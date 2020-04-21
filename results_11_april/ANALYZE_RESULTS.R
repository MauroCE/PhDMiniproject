library(tidyverse)

# SETTINGS
pop = "stable"
path = paste("~/Documents/University/PhDMiniproject/results_11_april/job_", pop, "_results/", sep="")


# GET TRUE VALUES
true_values <- c(10, 1, 1)
if (pop == "growing")     true_values <- c(20, 20, 0.05)
if (pop == "contracting") true_values <- c(1, 0.05, 1)

# FUNCTION TO READ IN A FILE
read_results <- function(file_name){
  unname(as.matrix(read.table(paste(path, file_name, sep=""))))
}

# FILES WITH TRACES FOR BOTH WORKERS (FULL AND THINNED TRACES)
w1_trace_full <- paste("x_history_full_w1_", pop, sep="")
w2_trace_full <- paste("x_history_full_w2_", pop, sep="")
w1_trace      <- paste("x_history_w1_", pop, sep="")
w2_trace      <- paste("x_history_w2_", pop, sep="")
# READ FILES
w1t  <- data.frame(read_results(w1_trace))
w2t  <- data.frame(read_results(w2_trace))
w1tf <- data.frame(read_results(w1_trace_full))
w2tf <- data.frame(read_results(w2_trace_full))
# REMOVE LAST ROW OF ZEROS
if (all(w1t[nrow(w1t), ] == rep(0, 3))) w1t <- w1t[1:(nrow(w1t)-1), ]
if (all(w2t[nrow(w2t), ] == rep(0, 3))) w2t <- w2t[1:(nrow(w2t)-1), ]
w1tf <- w1tf[1:(nrow(w1tf)-20), ]
w2tf <- w2tf[1:(nrow(w2tf)-20), ]
# PLOT TRACES OF ONE WORKER
plot_traces <- function(wtf){
  wtf %>% 
    mutate(rn=row_number()) %>% 
    gather("key", "value", -rn) %>% 
    {
      ggplot(data=. , aes(x=rn, y=value)) + 
        geom_line() + 
        facet_wrap(~key, ncol = 3)
    }
}
# FUNCTION TO PLOT ALL TRACES TOGETHER
plot_traces_together <- function(df1, df2){
  df1 <- mutate(df1, worker=as.factor(1), iteration=row_number())
  df2 <- mutate(df2, worker=as.factor(2), iteration=row_number())
  df  <- rbind(df1, df2)
  df <- gather(df, "key", "trace", -worker, -iteration)
  ggplot(data=df) + geom_line(aes(x=iteration, y=trace, color=worker)) + facet_wrap(~key) + ggtitle("Workers Traces") + 
    theme(plot.title = element_text(hjust=0.5))
}
# PLOT THEM
plot_traces_together(w1t, w2t)
plot_traces_together(w1tf, w2tf)
# THETA POSTERIOR HISTORY
tp_history     <- paste("tp_history_", pop, sep="")
tp_history_avg <- paste("tp_history_avg_", pop, sep="")
polyak_history <- paste("polyak_history_", pop, sep="")
# READ THEM
tp     <- unname(as.matrix(read_results(tp_history)))
tp_avg <- unname(as.matrix(read_results(tp_history_avg)))
polyak <- unname(as.matrix(read_results(polyak_history)))
# REMOVE ZEROS
tp     <- tp[1:(nrow(tp)-1), ]
tp_avg <- tp_avg[1:(nrow(tp_avg)-1), ]
polyak <- polyak[1:(nrow(polyak)-1), ]


# PLOT THEM
natural2muSigma <- function(theta, d){
  diagonal <- -1/theta[(d+1):(2*d)]
  return(list(mu=theta[1:d]*diagonal, sigmas=diagonal))
}
path_plot <- function(path_values, true_values){
  d <- 3
  # Remove last row if it contains only zeros
  n <- nrow(path_values)
  if (all(path_values[n, ] == rep(0.0, ))){
    path_values <- path_values[1:(n-1), ]
    n          <- n-1
  }
  # Find means and standard deviations
  means  <- matrix(0.0, nrow=n, ncol=d)
  sds    <- matrix(0.0, nrow=n, ncol=d)
  for (j in 1:n){
    mean_params <- natural2muSigma(matrix(path_values[j, ]), d)
    means[j, ] <- exp(mean_params$mu + 0.5*mean_params$sigmas)
    sds[j, ]   <- sqrt(diag(diag(means[j, ]) %*% (exp(diag(mean_params$sigmas)) - 1) %*% diag(means[j, ])))
  }
  # Compute data frames for plotting
  min_means <- means - 1.96 * sds
  max_means <- means + 1.96 * sds
  means_df <- means %>% data.frame()     %>% mutate(index=row_number()) %>% gather("Dimension", "Path", -index) 
  min_df   <- min_means %>% data.frame() %>% mutate(index=row_number()) %>% gather("Dimension", "Path", -index)
  max_df   <- max_means %>% data.frame() %>% mutate(index=row_number()) %>% gather("Dimension", "Path", -index)
  final_df <- mutate(means_df, ymin=min_df$Path, ymax=max_df$Path)
  cutoff <- data.frame(values=true_values, `Dimension`=c("X1", "X2", "X3"), stringsAsFactors = FALSE)
  # Plot 
  ggplot(data=final_df) + 
    geom_ribbon(aes(x=index, ymin=ymin, ymax=ymax), fill="grey70") + 
    geom_line(aes(x=index, y=Path)) + 
    facet_wrap(~ `Dimension`, nrow=1) +
    geom_hline(data=cutoff, aes(yintercept=values), lty=2, alpha=0.5) +
    ggtitle("Path of Each Coordinate") +
    theme(plot.title=element_text(hjust=0.5))
}

path_plot(tp, true_values)
path_plot(tp_avg, true_values)
path_plot(polyak, true_values) 

