# convert excel seasonly adjusted by eviews in R and take logs

library(haven)
library(readr)
library(readxl)
library(dplyr)
library(zoo)
library(stringr)

all <- read_excel("../data/adjusted_data_x.xlsx", sheet = 1, skip=0)

# get time in years
str(all)
all$time_y <- as.numeric(as.yearmon(all$time, format = "%Y %b"))

# remove SA in colnames
remove_sa <- str_detect(colnames(all),"_SA$") # find col nums with "_SA" in the end

old_names <- colnames(all)[remove_sa]
new_names <- str_match(old_names,"(.+)_SA")[,2]
colnames(all)[remove_sa] <- new_names

colnames(all)

write_csv(all, "../data/data_2015_after_eviews.csv")


df_sa <- read_csv("../data/data_2015_after_eviews.csv")
colnames(df_sa)

df_final <- mutate_each(df_sa, "log", -ib_rate, -lend_rate, -unemp_rate, -gov_balance, -time_y)
head(df_final)
head(df_sa)

write_csv(df_final, "../data/df_2015_eviews_ln.csv")


