
library(tidyverse)
library(reshape2)
library(rio)
library(stringr)


old_table <- import('old_table_replica.odf.ods')
colnames(old_table)[1] <- 'variable'
old_table_model <- import('old_table_replica_model.ods')
colnames(old_table_model)[1] <- 'variable'


rel_rmse_table <- read_rds('rel_rms_table.Rds')



# univariate leader  ------------------------------------------------------
univariate <- rel_rmse_table %>% rename(h = row_number) %>%
  filter(model_type %in% c('arima', 'ets', 'rw')) %>% select(-rmse, -mae) %>% filter(var_set == 'set_C') %>%
  select(-var_set, -pars_id) %>% filter(h %in% c(1, 3, 6, 9, 12))

uni_leader_long <- univariate %>% group_by(h, variable) %>% top_n(-1, mse) %>% ungroup() 


uni_leader_long %>% dcast(variable ~ h, value.var = 'model_type')
uni_leader_long %>% dcast(variable ~ h, value.var = 'mse')


ggplot(data = uni_leader_long, aes(x = factor(h), y = variable)) +
  geom_tile(aes(fill = model_type)) +
  geom_text(aes(label = round(mse, 2))) +
  theme_minimal() + labs(x = 'Forecasting horizon', 
                         title = 'Univariate models MSE relative to random walk MSE',
                         y = 'Variable', 
                         fill = 'Model')

# conclusion:
# random walk is a good benchmark model
# but not always the perfect one


# bad quality of var-lasso ------------------------------------------------
var_lasso <- rel_rmse_table %>% rename(h = row_number) %>%
  filter(model_type %in% c('var_lasso')) %>% select(-rmse, -mae) %>%
  select(-model_type, -pars_id) %>% filter(h %in% c(1, 3, 6, 9, 12))

var_lasso_leader_subset <- var_lasso %>% group_by(h, variable) %>% top_n(-1, mse) %>% ungroup() 

ggplot(data = var_lasso_leader_subset, aes(x = factor(h), y = variable)) +
  geom_tile(aes(fill = var_set)) +
  geom_text(aes(label = round(mse, 2))) +
  theme_minimal() + labs(x = 'Forecasting horizon', 
                         title = 'VAR-LASSO MSE relative to random walk MSE',
                         y = 'Variable', 
                         fill = 'Variable set')

# or we may just say, that the quality of var-lasso is bad (hide these 2000+ values)
# impovement for only one variable: gov_balance
# as comment: almost always set_C is better


# var-lasso for gov balance -----------------------------------------------
# unique series where var-lasso is the best :)
all_gov_balance <- rel_rmse_table %>% rename(h = row_number) %>% filter(variable == 'gov_balance') %>%
  filter(model_type %in% c('var_lasso', 'arima', 'ets')) %>% select(-rmse, -mae) %>%
 filter(h %in% c(1, 3, 6, 9, 12), var_set == 'set_C') 

all_gov_4_plot <- all_gov_balance %>% select(-variable, -var_set) %>% group_by(h, model_type) %>%
  top_n(-1, mse) %>% ungroup()

ggplot(all_gov_4_plot, aes(x = h, y = mse)) + geom_line(aes(color = model_type)) +
  labs(x = 'Forecasting horizon',  title = 'Model MSE relative to random walk MSE for gov_balance',
  y = 'MSE relative to random walk',  color = 'Model type')


# final leader ------------------------------------------------------------
old_table_melt <- melt(old_table, id.vars = 'variable', value.name = 'mse', variable.name = 'h') %>%
   mutate(h = as.numeric(str_sub(h, start = 2)))
old_table_melt


old_table_mod_melt <- melt(old_table_model, id.vars = 'variable', value.name = 'model_type', variable.name = 'h') %>%
  mutate(h = as.numeric(str_sub(h, start = 2)), model_type = str_to_lower(model_type))
old_table_mod_melt

old_table_restored <- left_join(old_table_melt, old_table_mod_melt, by = c('variable', 'h'))
old_table_restored

new_leader <- rel_rmse_table %>% rename(h = row_number) %>%
  filter(model_type %in% c('arima', 'ets', 'rw', 'var_lasso')) %>% select(-rmse, -mae) %>% filter(var_set == 'set_C') %>%
  select(-var_set) %>% filter(h %in% c(1, 3, 6, 9, 12)) %>%
  group_by(h, variable, model_type) %>% top_n(-1, mse) %>% ungroup() %>% select(-pars_id)
new_leader


all_models <- bind_rows(old_table_restored, new_leader)
all_models_lead <- group_by(all_models, h, variable) %>% top_n(-1, mse) %>% ungroup()


all_models_lead 
ggplot(data = all_models_lead, aes(x = factor(h), y = variable)) +
  geom_tile(aes(fill = model_type)) +
  geom_text(aes(label = round(mse, 2))) +
  theme_minimal() + labs(x = 'Forecasting horizon', 
                         title = 'Model MSE relative to random walk MSE',
                         y = 'Variable', 
                         fill = 'Model')
