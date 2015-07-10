# test 1 model estimation

source("400_model_funs.R")

mlist <- read_csv("../estimation/mlist_A.csv")
model_info <- mlist %>% filter(id==3)

estimate_model(model_info)
