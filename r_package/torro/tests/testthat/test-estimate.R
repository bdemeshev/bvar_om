context("Estimating models")

test_that("estimate_rw dummy test", {
  y_small <- rus_macro[, c("cpi", "employment")]
  rw_model <- estimate_rw(y_small)
  expect_equal(class(rw_model), "character")
})


test_that("estimate_arima test", {
  y_small <- rus_macro[, c("cpi", "employment")]
  arima_model <- estimate_arima(y_small)
  expect_equal(class(arima_model), "list")
  expect_equal(class(arima_model[[1]]), c("ARIMA", "Arima"))
})

test_that("estimate_var test", {
  y_small <- rus_macro[, c("cpi", "employment")]
  var_model <- estimate_var(y_small)
  expect_equal(class(var_model), "varest")
})