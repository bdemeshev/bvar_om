# small parallel demo

# it's a good idea to use plain R (not Rstudio) for parallel computations

# this function will be used by all cores:
my_fun <- function(x) {
  a <- ggplot2::qplot(rnorm(100))
  return(x^2)
}

library(doParallel)  


export_functions <- c("my_fun")
export_packages <- c("ggplot2", "forecast", "BigVAR")

n_clusters <- 4
# one may try to detect cores
detectCores()
# problems:
# 1. this function is very platform dependent and may give wrong results 
# better works for mac and linux
# 2. optimal number of cores is usually lower than the total number of cores
# so the optimal number of cores may be determined by experiments

# "cat" and "message" function will be redirected to "log.txt"
cluster <- makeCluster(n_clusters, outfile = "log.txt")

registerDoParallel(cluster)

foreach(i = 1:10000, 
        .export = export_functions, 
        .packages = export_packages) %dopar% {
  r <- my_fun(i)
}

stopCluster(cluster)  # a must
