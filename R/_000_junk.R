library(BigVAR)

data(Y)
y_st <- apply(Y, 2, FUN = scale)
y_st_ts <- ts(y_st, start = 1, frequency = 1)

model_spec <- constructModel(y_st, p = 1, 
                             struct = "OwnOther", gran = c(10, 2))
model_est <- BigVAR.est(model_spec)

model_cv <- cv.BigVAR(model_spec)


devtools::install_github("wbnicholson/BigVAR/BigVAR")


library(bvarsv)
model <- bvar.sv.tvp(y_subset, p = 2, nrep = 5000, nburn = 1000)
mean(predictive.draws(model, v = 2, h = 2)$y)
