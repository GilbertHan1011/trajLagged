test_that("sumWeight works correctly", {
  # Test basic functionality
  x <- c(1, 2, 3, 4, 5)
  weight <- c(0.1, 0.2, 0.3, 0.2, 0.2)
  
  result <- sumWeight(x, weight)
  
  # Should return a numeric value
  expect_type(result, "double")
  expect_length(result, 1)
  
  # Should handle edge cases
  expect_no_error(sumWeight(c(1), c(1)))
  expect_no_error(sumWeight(c(1, 1, 1), c(1, 1, 1)))
})

test_that("get_value_gam works correctly", {
  # Create mock data frame
  df_gam <- data.frame(
    time = seq(0, 1, length.out = 10),
    reference_smoothed = seq(1, 2, length.out = 10),
    target_smoothed = seq(1.5, 2.5, length.out = 10),
    sum_weight = seq(1, 0, length.out = 10)
  )
  
  result <- get_value_gam(df_gam)
  
  # Should return a numeric value
  expect_type(result, "double")
  expect_length(result, 1)
  
  # Target should be different from reference
  expect_true(is.finite(result))
}) 