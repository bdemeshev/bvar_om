library("foreach")

parallel <- "off" # "windows"/"unix"/"off"

source("400_model_funs.R")

# create model list
mlist <- create_model_list()
write_csv(mlist, path = "../estimation/mlist_A.csv")

# estimate models from list
mlist <- read_csv("../estimation/mlist_A.csv")
mlist <- estimate_models(mlist, parallel = parallel) # status and filename are updated
write_csv(mlist, path = "../estimation/mlist_A.csv")

# create prediction list
# plist <- create_prediction_list(...)
# write_csv(plist, path = "../estimation/plist_A.csv")

# make predictions
# plist <- read_csv("../estimation/plist_A.csv")
# plist <- make_predictions(plist, parallel = parallel) # status and filename are updated
# write_csv(plist, path = "../estimation/plist_A.csv")

# evaluate models
# plist <- read_csv("../estimation/plist_A.csv")
# evlist <- evaluate_models(plist)



