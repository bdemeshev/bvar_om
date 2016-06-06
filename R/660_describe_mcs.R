# 660 describe MCS results

library("dplyr")
library("ggplot2")
library("MCS")
library("reshape2")
library("tidyr")
library("scales")
library("xtable")


base_set_A <- c("ind_prod", "cpi", "ib_rate")
base_set_B <- c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price")
base_set_C <- c("employment", "ind_prod", "cpi", 
                "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
                "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                "agriculture", "retail", "gov_balance", "export", "import")

all_h <- 1:12
all_variables <- c("employment", "ind_prod", "cpi", 
                   "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
                   "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                   "agriculture", "retail", "gov_balance", "export", "import")

h <- 2
analysed_variable <- "ind_prod"

mcs_folder <- "../estimation/tables_rmsfe/"



# b_x_y_z - BVAR
# x - var_set:
# 1 - set_A
# 2 - set_B
# 3 - set_C
# y - fit_set
# 1 - ind+cpi
# 2 - ind+cpi+rate
# z - number of lags in BVAR model

# v_x_y_z: - VAR/RW/WN
# v_1 - rw 
# v_2 - var
# v_3 - wn
# y - var_set:
# 1 - set_A
# 2 - set_B
# 3 - set_C
# for RW and WN var_set is set to set_C
# z - number of lags in a VAR model 

# for fixed h and analysed variable the number of models is 
# bvar: 72 = 2 fit_set * 3 var_set * 12 lags
# var: 26 = 1 wn + 1 rw + 2 var_set * 12 lags

# in comparison we choose 1 fit_set, so
# 36 bvars + 24 vars + 1 rw + 1 wn = 62 models are compared


all_survived <- NULL
for (analysed_variable in all_variables) {
  for (h in all_h) {
    filename <- paste0(mcs_folder, "best_models_", analysed_variable, "_", h, ".Rds" )

    mcs <- readRDS(filename)
    survived <- data.frame(model = rownames(mcs@show), mcs@show) 
    survived <- survived %>% arrange(Rank_M)
    survived <- mutate(survived, h = h, variable = analysed_variable)
    all_survived <- rbind(all_survived, survived)
  }
}

surv <- all_survived %>% filter(h == 3, variable == "ind_prod")
surv

random_walks <- all_survived %>% filter(model == "v_1_3_1")
random_walks

rw_table <- dcast(data = random_walks, variable ~ h, fun.aggregate = length)
rw_table
 
white_noises <- all_survived %>% filter(model == "v_3_3_0")
wn_table <- dcast(data = white_noises, variable ~ h, fun.aggregate = length)
wn_table

all_survived_sep <- all_survived %>%
  separate(model, into = c("model_type", "code1", "code2", "lag"), sep = "_")

table_b1 <- filter(all_survived_sep, model_type == "b", code1 == 1)
table_b1_grouped <- group_by(table_b1, h, variable) %>% summarise(n_survived = n())
table_b1_grouped <- mutate(table_b1_grouped, name = "b1")
b1_wide <- dcast(table_b1_grouped, variable ~ h, fill = 0, value.var = "n_survived")

table_b2 <- filter(all_survived_sep, model_type == "b", code1 == 2)
table_b2_grouped <- group_by(table_b2, h, variable) %>% summarise(n_survived = n())
table_b2_grouped <- mutate(table_b2_grouped, name = "b2")
b2_wide <- dcast(table_b2_grouped, variable ~ h, fill = 0, value.var = "n_survived")

table_b3 <- filter(all_survived_sep, model_type == "b", code1 == 3)
table_b3_grouped <- group_by(table_b3, h, variable) %>% summarise(n_survived = n())
table_b3_grouped <- mutate(table_b3_grouped, name = "b3")
b3_wide <- dcast(table_b3_grouped, variable ~ h, fill = 0, value.var = "n_survived")


table_v1 <- filter(all_survived_sep, model_type == "v", code1 == 2, code2 == 1)
table_v1_grouped <- group_by(table_v1, h, variable) %>% summarise(n_survived = n())
table_v1_grouped <- mutate(table_v1_grouped, name = "v1")
v1_wide <- dcast(table_v1_grouped, variable ~ h, fill = 0, value.var = "n_survived")

table_v2 <- filter(all_survived_sep, model_type == "v", code1 == 2, code2 == 2)
table_v2_grouped <- group_by(table_v2, h, variable) %>% summarise(n_survived = n())
table_v2_grouped <- mutate(table_v2_grouped, name = "v2")
v2_wide <- dcast(table_v2_grouped, variable ~ h, fill = 0, value.var = "n_survived")

all_grouped <- bind_rows(table_b1_grouped, table_b2_grouped, table_b3_grouped, table_v1_grouped, table_v2_grouped)
all_wide <- dcast(all_grouped, variable ~ h + name, fill = 0, value.var = "n_survived")

margin_variable <- group_by(all_grouped, variable, name) %>% summarise(ratio_survived =  round(sum(n_survived)/(12*12),2))
margin_var_wide <- dcast (margin_variable, variable ~ name, value.var = "ratio_survived")

margin_h <- group_by(all_grouped, h, name) %>% summarise(ratio_survived =  round(sum(n_survived)/(23*12),2))
margin_h_wide <- dcast(margin_h,h ~ name, value.var = "ratio_survived")
margin_h_wide
xtable(margin_h_wide) # transform to Latex

# some plots
one_selection <- all_grouped %>% filter(h == 5, variable == "ind_prod")

one_histo <- ggplot(data = one_selection) + 
  geom_bar(aes(x = name, y = n_survived), stat = "identity") +
  scale_y_continuous(breaks = seq(from = 2, to = 12, by = 2), 
                     labels = seq(from = 2, to = 12, by = 2)) 
one_histo


all_selection <- all_grouped %>% filter(h %in% c(1, 3, 6, 9, 12),
                                        variable %in% base_set_B)
all_histo <- ggplot(data = all_selection) + 
  geom_bar(aes(x = name, y = n_survived), stat = "identity") +
  scale_y_continuous(breaks = seq(from = 2, to = 12, by = 2), 
                     labels = seq(from = 2, to = 12, by = 2)) +
  facet_grid(variable ~ h)
all_histo



