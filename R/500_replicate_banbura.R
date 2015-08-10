# 500_replicate_banbura.R

# input: "../data/df_2015_final.csv", "../data/var_set_info.csv"
# output: 



library("foreach")

source("400_model_funs.R")

parallel <- "off" # "windows"/"unix"/"off"
ncpu <- 30

