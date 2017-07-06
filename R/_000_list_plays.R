# play with model/var_set/shifts lists

source("_000_forecast_mts.R")

rus_macro <- load_rus_data()

shifts <- tribble(~shift_name, ~T_start, ~win_type, ~win_start_length, ~n_shifts,
                  "base", 1, "moving", 120, 10)

create_var_set_info <- function() {
  add_A <- dplyr::tibble(var_set = "set_A", variable = c("ind_prod", "cpi", "ib_rate"))
  add_B <- dplyr::tibble(var_set = "set_B", variable = c("ind_prod", "cpi", "ib_rate", "m2", "reer", "oil_price"))
  add_C <- dplyr::tibble(var_set = "set_C", variable = c("employment", "ind_prod", "cpi", 
                                                             "ib_rate", "lend_rate", "real_income", "unemp_rate", "oil_price", "ppi", "construction", 
                                                             "real_investment", "wage", "m2", "reer", "gas_price", "nfa_cb", "ner", "labor_request", 
                                                             "agriculture", "retail", "gov_balance", "export", "import"))
  
  var_set_info <- dplyr::bind_rows(add_A, add_B, add_C)
  var_set_info <- dplyr::mutate(var_set_info, pre_transform = "none", post_transfomr = "none")
  
  return(var_set_info)
}

vars_sets <- create_var_set_info()



