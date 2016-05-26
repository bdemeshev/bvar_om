library("dplyr")
library("reshape2")
library("readr")

data_for_tables_file <- "../estimation/tables_rmsfe/tables_rmsfe_raw.Rds"
final_table_file <- "../estimation/tables_rmsfe/table_rmsfe_final.csv"


raw_t_data <- readRDS(data_for_tables_file)

raw2 <- filter(raw_t_data, variable == analysed, 
               fit_set == "ind+cpi+rate") %>% select(-analysed, -fit_set)

raw_melted <- melt(raw2, id.vars = c("h", "n_lag", "variable"), variable.name = "model")

best_model <- group_by(raw_melted, h, variable) %>% filter(value == min(value))


table_rmsfe <- dcast(best_model, variable ~ h, value.var = "value")
table_final <- select(table_rmsfe, variable, `1`, `3`, `6`, `9`, `12`)

table_rmsfe_m <- dcast(best_model, variable ~ h, value.var = "model")
table_final_m <- select(table_rmsfe_m, variable, `1`, `3`, `6`, `9`, `12`)

write_csv(table_final_m, file = final_table_file)

#### test for equaility of VAR and BVAR for set A
best_model %>% group_by(h, variable) %>% summarise(n = n()) %>% arrange(-n)
ag <- raw_melted %>% filter(variable == "agriculture", h == 1) %>% arrange(value)
head(ag)
ag$value[1] == ag$value[2]
# moral: var and bvar for set A are equal up to machine precision
