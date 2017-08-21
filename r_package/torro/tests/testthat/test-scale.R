context("Scaling time series")

test_that("get_scales just calculates mean and sd", {
  expect_equal(get_scales(cars)$mu[1], 42.98)
  expect_equal(get_scales(cars)$sd[1], 25.76937749)
})


test_that("scale_to original scale should do nothing", {
  cars_scaled <- scale_to(cars, get_scales(cars))
  expect_equal(cars_scaled$dist, cars$dist)
  expect_equal(cars_scaled$speed, cars$speed)
})