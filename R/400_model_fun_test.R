# test 1 model estimation

source("400_model_funs.R")

df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")
var_set_info <- var_set_info[-11,]

mlist <- read_csv("../estimation/mlist_A.csv")
model_info <- mlist %>% filter(id==3)
model_info
estimate_model(model_info)

# тестируем выбрасывание переменной
for (j in 10:31) {
  for (i in (j+1):32) {
    var_set_info <- read_csv("../data/var_set_info.csv")
    var_set_info <- var_set_info[-c(i,j),]
    estimate_model(model_info)
  }
}



