# 200_load_after_eviews

# convert excel seasonly adjusted by eviews in R and take logs

# input: "../data/adjusted_data_x.xlsx"
# output: "../data/data_2015_after_eviews.csv" (just correction of colnames)
# output: "../data/df_2015_final.csv" (filter from 1995 year, take logs)

library(haven)
library(readr)
library(readxl)
library(dplyr)
library(zoo)
library(stringr)

all <- read_excel("../data/adjusted_data_x.xlsx", sheet = 1, skip=0)

# get time in years
# str(all)
# all$time_y <- as.numeric(as.yearmon(all$time, format = "%Y %b"))

# remove SA in colnames
remove_sa <- str_detect(colnames(all),"_SA$") # find col nums with "_SA" in the end

old_names <- colnames(all)[remove_sa]
new_names <- str_match(old_names,"(.+)_SA")[,2]
colnames(all)[remove_sa] <- new_names

# remove 2 in colnames
remove_2 <- str_detect(colnames(all),"2$")
old_names <- colnames(all)[remove_2]
new_names <- str_match(old_names,"(.+)2")[,2]
colnames(all)[remove_2] <- new_names



colnames(all) <- tolower(colnames(all))

colnames(all)[colnames(all)=="m2_constprice"] <- "m2"

write_csv(all, "../data/data_2015_after_eviews.csv")


df_sa <- read_csv("../data/data_2015_after_eviews.csv")
colnames(df_sa)

dplyr::select(df_sa,obs,time,time_y) %>% head()

df_sa <- dplyr::select(df_sa, -time, -obs)

df_final <- mutate_each(df_sa, "log", -ib_rate, -lend_rate, -unemp_rate, -gov_balance, -time_y)
head(df_final)
head(df_sa)

df_final <- filter(df_final, time_y >= 1995)

write_csv(df_final, "../data/df_2015_final.csv")


