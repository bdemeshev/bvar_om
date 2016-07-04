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



### get survived random walks
full_table <- expand.grid(h = 1:12, variable = base_set_C, stringsAsFactors = FALSE)


random_walks <- all_survived %>% filter(model == "v_1_3_1")
random_walks$survived <- 1

dead_rw <- anti_join(full_table, random_walks, by = c("h", "variable"))
dead_rw$survived <- 0

random_walks <- bind_rows(random_walks, dead_rw)

rw_table <- dcast(data = random_walks, variable ~ h, value.var = "survived")
rw_table


### get survived white noises
full_table <- expand.grid(h = 1:12, variable = base_set_C, stringsAsFactors = FALSE)

white_noises <- all_survived %>% filter(model == "v_3_3_0")
white_noises$survived <- 1

dead_wn <- anti_join(full_table, white_noises, by = c("h", "variable"))
dead_wn$survived <- 0

white_noises <- bind_rows(white_noises, dead_wn)

wn_table <- dcast(data = white_noises, variable ~ h, value.var = "survived")
wn_table


### joint rw+wn table

wn_table_4join <- dplyr::select(white_noises, h, variable, wn_survived = survived)
rw_table_4join <- dplyr::select(random_walks, h, variable, rw_survived = survived)
rwwn_table <- dplyr::left_join(wn_table_4join, rw_table_4join, by = c("h", "variable"))

rwwn_table <- mutate(rwwn_table, 
                     wn_letter = ifelse(wn_survived == 1, "WN", ""),
                     rw_letter = ifelse(rw_survived == 1, "RW", ""),
                     label = paste0(rw_letter, wn_letter)
                     )


### percentage of survived models

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

### all histo at once
selected_h <- c(1, 3, 6, 9, 12)
selected_variables <- setdiff(base_set_B, "m2")

all_selection <- all_grouped %>% filter(h %in% selected_h,
                                        variable %in% selected_variables)

# remove small BVAR as it is equal to VAR:
all_selection <- filter(all_selection, !name == "b1")




model_transform <- function(x) {
  
  model_names <- list(
    'b2' = 'BVAR6/7',
    'b3' = 'BVAR23',
    'v1' = 'VAR3/4',
    'v2' = 'VAR6/7' 
  )
  
  return(unname(unlist(model_names[x])))
}


all_selection$name <- model_transform(all_selection$name)


h_labeller <- as_labeller(c('1' = 'h = 1',
                            '3' = 'h = 3',
                            '6' = 'h = 6',
                            '9' = 'h = 9',
                            '12' = 'h = 12'))

all_histo <- ggplot(data = all_selection) + 
  geom_bar(aes(x = name, y = n_survived), stat = "identity") +
  scale_y_continuous(breaks = seq(from = 2, to = 12, by = 2), 
                     labels = seq(from = 2, to = 12, by = 2)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(variable ~ h,  labeller = labeller(h = h_labeller)) + 
  ylab("Number of survived models") + xlab("")
all_histo 

# add rwwn tag in the right corner 
rwwn_table_selection <- rwwn_table %>% filter(h %in% selected_h,
                                              variable %in% selected_variables)
rwwn_table_selection <- mutate(rwwn_table_selection, x = Inf, y = Inf)
# x = Inf and y = Inf means right corner

all_histo_rwwn <- all_histo + 
              geom_text(aes(x, y, label = label), data = rwwn_table_selection, 
              hjust = 1, vjust = 1, size = 4)
all_histo_rwwn






