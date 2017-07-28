my_fun <- function(x) {
  a <- ggplot2::qplot(rnorm(100))
  return(x^2)
}

library(doParallel)  

export_functions <- c("my_fun")
export_packages <- c("ggplot2", "forecast", "BigVAR")

cluster <- makeCluster(4, outfile = "log.txt")
registerDoParallel(cluster)

foreach(i = 1:10000, 
        .export = export_functions, 
        .packages = export_packages) %dopar% {
  r <- my_fun(i)
}

stopCluster(cluster)  
