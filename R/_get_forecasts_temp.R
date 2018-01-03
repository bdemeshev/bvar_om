library(tidyverse)
library(torro)
library(stringr)

a <- list.files("~/Documents/bvar_om/estimation/bvar_alternatives/fits/")
df <- tibble(fnames = a)
df <- df %>% mutate(number = as.numeric(str_extract(fnames, "[0-9]+")))
missing <- setdiff(1:4995, df$number)

dput(missing)



ff_from <- "~/Documents/bvar_om/estimation/bvar_alternatives/fits/"

ff_to <- "~/Documents/bvar_om/estimation/refits/fits/"

files <- list.files(ff_from)

for (file in files) {
  cat(file, "\n")
  old_fit <- read_rds(paste0(ff_from, file))
  new_fit <- reforecast_fit(old_fit, 12)
  write_rds(new_fit, path = paste0(ff_to, file))
}


basefolder <- "~/Documents/bvar_om/estimation/refits/"
all_forecasts <- get_forecasts_from_fit_files(basefolder = basefolder)
write_rds(all_forecasts, paste0(basefolder, "all_forecasts.Rds"))
