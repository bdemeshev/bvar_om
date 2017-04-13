library(BigVAR)
library(forecast)
library(forecastHybrid)
library(caret)

# load data
df <- readr::read_csv("../data/df_2015_final.csv")
df <- dplyr::select(df, -time_y)

