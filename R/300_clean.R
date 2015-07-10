# deseason (clean data)
library(stringr)
library(readr)

user_name <- Sys.info()[7]

if (user_name=="boris") {
  Sys.setenv(X13_PATH ="~/Documents/.bin/")
}

if (str_detect(user_name,"xana")) {
  Sys.setenv(X13_PATH="D:/Oxana/Science/Projects/ProjectXIV/x13as/")
}


library(seasonal)
checkX13()

all <- read_csv("../data/data_2015.csv")

all_95 <- dplyr::filter(all, time_y >= 1995 ) 

# unselect trade_balance and rts
df <- all_95 %>% dplyr::select(-trade_balance,-rts)


