# create variable sets

# input: "../data/df_2015_sa.csv"
# output: "../data/var_set_info.csv"

library(readr)


df <- read_csv("../data/df_2015_final.csv")


colnames(df)
dput(colnames(df))

add_3 <- data_frame(var_set="set_3", variable=c("ind_prod", "cpi", "ib_rate"))
add_6 <- data_frame(var_set="set_6", variable=c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))

# set with all 23 variables
add_23 <- data_frame(var_set="set_23",variable=c("employment", "ind_prod", "cpi", "ib_rate", "lend_rate", "real_income", 
                                                 "unemp_rate", "oil_price", "ppi", "construction", "real_investment", 
                                                 "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                                                 "agriculture", "retail", "gov_balance", "export", "import") )

var_set_info <- rbind(add_3,add_6,add_23)

write_csv(var_set_info,"../data/var_set_info.csv")
