# convert excel and other junk to Rds (or csv)

library(haven)
library(readr)
library(readxl)
library(dplyr)
library(zoo)


all <- read_excel("../data/data_bvar_2015.xlsx", sheet = 1, skip=5)

# get time in years
all$time_y <- as.numeric(as.yearmon(all$time, format = "%Y %b"))

colnames(all)

write_csv(all, "../data/data_2015.csv")


