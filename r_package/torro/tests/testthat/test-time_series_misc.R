context("Time series misc functions")

test_that("next_obs_time gives first forecast time period", {
  y <- ts(1:5, start = c(2015, 2), frequency = 12)
  expect_equal(next_obs_time(y), 2015.5)
})


test_that("mforecast_to_matrix", {
  data(rus_macro)
  y_small <- rus_macro[, c("cpi", "employment", "m2")]
  var_model <- vars::VAR(y_small)
  y_for <- forecast::forecast(var_model, h = 2)
  ans <- matrix(c(5.00534477847833, 5.00716061746032, 4.55412742698148, 
                  4.55487075267483, 5.37511778301469, 5.38072113452549), nrow = 2)
  colnames(ans) <- c("cpi", "employment", "m2")
  expect_equal(mforecast_to_matrix(y_for), ans)
})


test_that("matrix to mforecast and back", {
  y_before <- matrix(1:6, nrow = 2)
  y_new <- matrix(11:16, nrow = 2)
  mfor <- matrix_to_mforecast(y_new, y_before, method = "junk")
  expect_equal(mforecast_to_matrix(mfor), y_new)
})