library(BigVAR)
library(forecast)
library(forecastHybrid)
library(caret)
library(Rtsne)
library(tidyverse)
library(ggrepel)


# load data
df <- readr::read_csv("../data/df_2015_final.csv")
df <- dplyr::select(df, -time_y)


df_scaled <- apply(df, MARGIN = 2, FUN = scale)

ts_data <- data.frame(t(df_scaled))
help(package = "Rtsne")

set.seed(420)
tsne_out <- Rtsne(ts_data, 
                  perplexity = 1, 
                  verbose = TRUE,
                  theta = 0.0,
                  max_iter = 10000)
# lower theta = higher accuracy, slower speed


plot_df <- data_frame(tsne_x = tsne_out$Y[, 1],
                      tsne_y = tsne_out$Y[, 2],
                      name = colnames(df))

qplot(data = plot_df, geom = "point", x = tsne_x, y = tsne_y) +
  geom_text_repel(aes(label = name))



# old good pca 

pca_out <- prcomp(ts_data, scale. = FALSE)
plot(pca_out, type = "l")

plot_df <- data_frame(name = colnames(df),
                      pca_1 = pca_out$x[, 1],
                      pca_2 = pca_out$x[, 2])


qplot(data = plot_df, geom = "point", x = pca_1, y = pca_2) +
  geom_text_repel(aes(label = name))






