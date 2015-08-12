
source("400_model_funs.R")

#### test 1: estimate big model (23 vars, 12 lags)


df <- read_csv("../data/df_2015_final.csv")
var_set_info <- read_csv("../data/var_set_info.csv")
#var_set_info <- var_set_info[-11,]

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_test.csv")

mlist <- read_csv("../estimation/mlist_test.csv")
model_info <- mlist %>% filter(id==3)
model_info
estimate_model(model_info)

#### test 2: estimate big model (23 vars, 12 lags) without pairs of variables

df <- read_csv("../data/df_2015_final.csv")

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_test.csv")

mlist <- read_csv("../estimation/mlist_test.csv")
model_info <- mlist %>% filter(id==3)

# тестируем выбрасывание переменной
for (j in 10:31) {
  for (i in (j+1):32) {
    var_set_info <- read_csv("../data/var_set_info.csv")
    var_set_info <- var_set_info[-c(i,j),]
    estimate_model(model_info)
  }
}

#### test 3: RW and WN prediction test

mlist <- create_rwwn_list()
plist <- data.frame(model_id=c(1,2), h=NA, type="in-sample")
mlist

pred_info <- plist[1,]
pred_info
forecast_model(pred_info, mlist, parallel = "off")

pred_info <- plist[2,]
pred_info
forecast_model(pred_info, mlist, parallel = "off")

