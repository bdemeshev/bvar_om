# create variable sets

# input: "../data/df_2015_sa.csv"
# output: "../data/var_set_info.csv"

library(readr)


df <- read_csv("../data/df_2015_sa.csv")


colnames(df)

add_3 <- data_frame(var_set="set_3", variable=c("ind_prod", "cpi", "ib_rate"))
add_6 <- data_frame(var_set="set_6", variable=c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))

# set with all 24 variables
add_24 <- data_frame(var_set="set_24",variable=colnames(df))

var_set_info <- rbind(add_3,add_6,add_24)

write_csv(var_set,"../data/var_set_info.csv")
