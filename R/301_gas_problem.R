# 301_gas_problem

library("readr")
library("zoo")
library("lubridate")
library("dplyr")

all <- read_csv("../data/data_2015.csv")

all_95 <- dplyr::filter(all, time_y >= 1995 ) 

# unselect trade_balance and rts
df <- all_95 %>% dplyr::select(-trade_balance, -rts)

gas_only <- df %>% select(time, gas_price)

glimpse(gas_only)
gas_only$time <- format(as.Date(as.yearmon(gas_only$time, format = "%Y %b")), "%Y-%m")

write_csv(gas_only, "../data/gas_only.csv", col_names = FALSE)

# upload to 
# http://www.seasonal.website/

