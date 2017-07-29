# small parallel demo

# this function will be used by all cores:
my_fun <- function(x) {
  a <- ggplot2::qplot(rnorm(100))
  return(x^2)
}

library(doParallel)  


export_functions <- c("my_fun")
export_packages <- c("ggplot2", "forecast", "BigVAR")

n_clusters <- 4

# "cat" and "message" function will be redirected to "log.txt"
cluster <- makeCluster(n_clusters, outfile = "log.txt")

registerDoParallel(cluster)

foreach(i = 1:10000, 
        .export = export_functions, 
        .packages = export_packages) %dopar% {
  r <- my_fun(i)
}

stopCluster(cluster)  
