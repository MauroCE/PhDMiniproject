library(tidyverse)

# SETTINGS
pop = "growing"

# FUNCTION TO READ IN A FILE
read_results <- function(file_name, path="~/Documents/University/PhDMiniproject/job_growing_results/"){
  unname(as.matrix(read.table(paste(path, file_name, sep=""))))
}

# GRAB RESULTS
file_name <- "x_history_full_w1_growing"
file_results <- unname(as.matrix(read.table(paste("~/Documents/University/PhDMiniproject/job_growing_results/", file_name, sep=""))))


# TRACES FOR WORKERS
w1_trace_full <- "x_history_full_w1_growing"
w2_trace_full <- "x_history_full_w2_growing"
w1_trace      <- "x_history_w1_growing"
w2_trace      <- "x_history_w2_growing"
# TODO: NEED TO REMOVE THE LAST ROWS OF ZEROS
w1t  <- data.frame(read_results(w1_trace))
w2t  <- data.frame(read_results(w2_trace))
w1tf <- data.frame(read_results(w1_trace_full))
w2tf <- data.frame(read_results(w2_trace_full))

# FUNCTION TO PLOT ALL TRACES OF ONLY ONE WORKER
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

# FUNCTION TO PLOT TRACES OF ALL WORKERS, COMPARED

