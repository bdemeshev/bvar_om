# test 1 model estimation

source("400_model_funs.R")

df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")

mlist <- read_csv("../estimation/mlist_A.csv")
model_info <- mlist %>% filter(id==3)

estimate_model(model_info)
