# model selection set

library("dplyr")
library("ggplot2")
library("MCS")

working_folder <- "../estimation/tables_rmsfe/"
forecast_filename_prefix <- "forecasts_"

all_vars <- c("employment", "ind_prod", "cpi", 
                "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
                "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                "agriculture", "retail", "gov_balance", "export", "import")

analysed_variable <- all_vars[2]
analysed_variable

filename <- paste0(working_folder, forecast_filename_prefix, analysed_variable, ".Rds")


d <- readRDS(filename)

str(d)

# d$actual_obs$variable

T_min <- 133
T_max <- max(d$actual_obs$t)

actual <- d$actual_obs %>% filter(t %in% T_min:T_max, 
                                  variable == analysed_variable)
bvar <- d$bvar_out_forecasts %>% filter(t %in% T_min:T_max, 
                                        variable == analysed_variable)
var <- d$rwwn_var_out_forecasts %>% filter(t %in% T_min:T_max, 
                                           variable == analysed_variable)

t_count_var <- var %>% group_by(t) %>% summarize(count = n())
t_count_var %>% head(15)
t_count_bvar <- bvar %>% group_by(t) %>% summarize(count = n())
t_count_bvar %>% head(15)

# total number of models (for t >= 144)
# bvar: 864
# var + rwwn: 312
# we have h=1, ..., 12
# for fixed h the number of models is 
# bvar: 72
# var: 26








