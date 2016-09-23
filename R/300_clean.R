# 300_clean
# deseason (clean data)
library(stringr)
library(readr)
library(lattice)

# input: "../data/data_2015.csv"
# output: "../data/df_2015_sa.csv", "../data/df_2015_final.csv"

user_name <- Sys.info()[7]

if (user_name == "boris") {
  Sys.setenv(X13_PATH = "~/Documents/.bin/")
}

if (str_detect(user_name, "xana")) {
  Sys.setenv(X13_PATH = "D:/Oxana/Science/Projects/ProjectXIV/x13as/")
}


library(seasonal)
checkX13()



seas_adjust <- function(df_ts, to_adjust = colnames(df_ts), 
                        method=c("X13", "stl"))  {
  method <- match.arg(method)
  df_sa <- df_ts # make a copy
  for (i in 1:ncol(df_ts)) {
    if (colnames(df_ts)[i] %in% to_adjust) {
      # cat(names(df_notime)[i], i, "\n") 
      message("Adjusting ", colnames(df_ts)[i])
      if (method == "X13") {
        model <- seas(df_ts[, i])
        df_sa[,i] <- final(model)
      }
      if (method == "stl") {
        res <- stl(df_ts[, i], s.window = "periodic")
        df_sa[, i] <- res$time.series[, 2] + res$time.series[, 3]
      } 
    }
  }
  return(df_sa)
}




all <- read_csv("../data/data_2015.csv")

all_95 <- dplyr::filter(all, time_y >= 1995 ) 

# unselect trade_balance and rts
df <- all_95 %>% dplyr::select(-trade_balance, -rts)

df_notime <- dplyr::select(df, -time)

psych::describe(df_notime)[, c("n", "mean","median", "sd")]


df_ts <- ts(df_notime, start = c(1995, 1), frequency = 12)
xyplot(df_ts)
# dput(colnames(df_ts))
# list of series to seasonally adjust
to_adjust <- c("employment", "ind_prod", "cpi", "real_income", 
               "unemp_rate", 
               "oil_price",  # maybe no?
               "ppi", "construction", "real_investment", 
               "wage", 
               "reer", # maybe no?
               # "gas_price", # ERROR adjusting gas_price
              "labor_request", 
               "agriculture", "retail", 
              "gov_balance", "export", "import"
) 

#to_adjust <- "gas_price"
df_sa <- seas_adjust(df_ts, to_adjust, method = "X13")

xyplot(df_sa)

write_csv(as.data.frame(df_sa),"../data/df_2015_sa.csv")



## take logs

df_sa <- read_csv("../data/df_2015_sa.csv")
colnames(df_sa)

df_final <- mutate_each(df_sa, "log", -ib_rate, -lend_rate, -unemp_rate, -gov_balance, -time_y)
head(df_final)
head(df_sa)

write_csv(df_final, "../data/df_2015_final.csv")







