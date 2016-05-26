# 500_banbura_usa_data.R

# input: "../data/df_2015_final.csv", "../data/var_set_info.csv"
# output: 

source("400_model_funs.R")
source("400_model_lists.R")
source("500_banbura_funs.R")
source("500_replicate_banbura_function.R")

# need to run only once
# source("200_load_after_eviews.R")


df <- read_csv("../data/usa_carriero.csv")
# названия переменных неизвестны:

add_3 <- data_frame(var_set = "set_3", variable = c("var1", "var2", "var3"))
add_6 <- data_frame(var_set = "set_6", variable = c("var1", "var2", "var3", "var4", "var5", "var6"))
add_13 <- data_frame(var_set = "set_13", variable = c("var1", 
                                                 "var2", 
                                                 "var3", 
                                                 "var4", "var5", "var6", 
                                                 "var7", "var8", 
                                                 "var9", 
                                                 "var10", "var11", 
                                                 "var12", 
                                                 "var13") )

var_set_info <- bind_rows(add_3,add_6,add_13)

fit_set_2vars <- data_frame(variable = c("var1", "var2"), fit_set = "var1+var2")
fit_set_3vars <- data_frame(variable = c("var1", "var2", "var3"), fit_set = "var1+var2+var3")

fit_set_info <- bind_rows(fit_set_2vars, fit_set_3vars)
fit_set_info

# error!
usa_data <- replicate_banbura(df, var_set_info = var_set_info, target_var_set = "set_3",
                              fit_set_info = fit_set_info) # need to check!!!
# причина ошибки: все названия var_set записаны жестко в 400_model_lists
# надобно их гибко создавать ;)


