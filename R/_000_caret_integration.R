library(BigVAR)
library(forecast)
library(forecastHybrid)
library(caret)

# load data
df <- readr::read_csv("../data/df_2015_final.csv")
df <- dplyr::select(df, -time_y)


p <- 12 # number of lags
window_width_dependent <- 120 # number of observations used for in-sample forecast in each window
window_width <- p + window_width_dependent # total window width (fixed)

h <- 6 # forecasting horizon

T_all <- nrow(df) # all available obs

TimeControl <- trainControl(method = "timeslice",
                              initialWindow = window_width,
                              horizon = h,
                              fixedWindow = TRUE)
str(TimeControl) # хрень какая-то


# sliced time indexes
time_slices <- createTimeSlices(1:T_all, 
                  initialWindow = window_width, 
                  horizon = h, 
                  fixedWindow = TRUE)

time_slices$train[[1]]
time_slices$test[[1]]

bigvar_caret <- list(type = "regression", loop = FALSE, library = "BigVAR")
prm <- dplyr::data_frame(parameter = c("lambda"),
                  class = c("numeric"),
                  label = c("Sigma (penalty)"))
# caret для одномерных только?
# увы!






