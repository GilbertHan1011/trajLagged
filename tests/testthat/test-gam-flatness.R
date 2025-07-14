test_that("test_gam_flatness works with flat data", {
  skip_if_not_installed("mgcv")
  
  # Create flat test data
  set.seed(123)
  n <- 50
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.1)  # Flat with small noise
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test flatness
  result <- test_gam_flatness(gam_model, method = "summary")
  
  # Check structure
  expect_type(result, "list")
  expect_true("summary_test" %in% names(result))
  expect_true("overall_interpretation" %in% names(result))
  
  # Check that summary test returns expected structure
  summary_result <- result$summary_test
  expect_true("method" %in% names(summary_result))
  expect_true("p_values" %in% names(summary_result))
  expect_true("is_flat" %in% names(summary_result))
  expect_true("interpretation" %in% names(summary_result))
  
  # For flat data, we expect high p-values (though not guaranteed due to randomness)
  expect_type(summary_result$p_values, "double")
  expect_type(summary_result$is_flat, "logical")
})

test_that("test_gam_flatness works with non-flat data", {
  skip_if_not_installed("mgcv")
  
  # Create non-flat test data
  set.seed(123)
  n <- 100
  x <- seq(0, 2*pi, length.out = n)
  y <- sin(x) + rnorm(n, 0, 0.1)  # Clear sinusoidal pattern
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test flatness
  result <- test_gam_flatness(gam_model, method = "summary")
  
  # Check structure
  expect_type(result, "list")
  expect_true("summary_test" %in% names(result))
  
  # For clearly non-flat data, we expect low p-values
  summary_result <- result$summary_test
  expect_type(summary_result$p_values, "double")
  expect_type(summary_result$is_flat, "logical")
  
  # With strong signal, should detect non-flatness
  expect_true(summary_result$p_values < 0.05)
  expect_false(summary_result$is_flat)
})

test_that("test_gam_flatness works with all methods", {
  skip_if_not_installed("mgcv")
  
  # Create test data
  set.seed(123)
  n <- 50
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.2)
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test with all methods
  result <- test_gam_flatness(gam_model, method = "all")
  
  # Check that all expected methods are present
  expected_methods <- c("summary_test", "wald_contrasts", "likelihood_ratio", "f_test")
  for (method in expected_methods) {
    expect_true(method %in% names(result), 
                info = paste("Method", method, "not found in results"))
  }
  
  # Check overall interpretation
  expect_true("overall_interpretation" %in% names(result))
  expect_type(result$overall_interpretation, "character")
})

test_that("test_gam_flatness handles individual methods correctly", {
  skip_if_not_installed("mgcv")
  
  # Create test data
  set.seed(123)
  n <- 50
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.1)
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test each method individually
  methods <- c("summary", "wald_contrasts", "likelihood_ratio", "f_test")
  
  for (method in methods) {
    result <- test_gam_flatness(gam_model, method = method)
    
    # Should contain the specific method and overall interpretation
    expect_true(paste0(method, "_test") %in% names(result) || method %in% names(result),
                info = paste("Method", method, "not found in results"))
    expect_true("overall_interpretation" %in% names(result))
  }
})

test_that("test_gam_flatness handles edge cases", {
  skip_if_not_installed("mgcv")
  
  # Test with invalid input
  expect_error(test_gam_flatness("not_a_gam"), 
               "gam_model must be a fitted GAM object")
  
  # Test with invalid method
  set.seed(123)
  n <- 20
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.1)
  data <- data.frame(x = x, y = y)
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  expect_error(test_gam_flatness(gam_model, method = "invalid_method"),
               "method must be one of")
})

test_that("Wald contrasts test works correctly", {
  skip_if_not_installed("mgcv")
  
  # Create test data
  set.seed(123)
  n <- 50
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.1)
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test Wald contrasts
  result <- test_gam_wald_contrasts(gam_model, n_points = 20)
  
  # Check structure
  expect_type(result, "list")
  expect_true("method" %in% names(result))
  expect_true("wald_statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("is_flat" %in% names(result))
  
  # Check types
  if (!is.na(result$wald_statistic)) {
    expect_type(result$wald_statistic, "double")
  }
  if (!is.na(result$p_value)) {
    expect_type(result$p_value, "double")
    expect_true(result$p_value >= 0 && result$p_value <= 1)
  }
})

test_that("Likelihood ratio test works correctly", {
  skip_if_not_installed("mgcv")
  
  # Create test data
  set.seed(123)
  n <- 50
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.1)
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test likelihood ratio
  result <- test_gam_likelihood_ratio(gam_model)
  
  # Check structure
  expect_type(result, "list")
  expect_true("method" %in% names(result))
  expect_true("lr_statistic" %in% names(result))
  expect_true("p_value" %in% names(result))
  expect_true("is_flat" %in% names(result))
  
  # Check types
  if (!is.na(result$lr_statistic)) {
    expect_type(result$lr_statistic, "double")
    expect_true(result$lr_statistic >= 0)
  }
  if (!is.na(result$p_value)) {
    expect_type(result$p_value, "double")
    expect_true(result$p_value >= 0 && result$p_value <= 1)
  }
})

test_that("F-test works correctly", {
  skip_if_not_installed("mgcv")
  
  # Create test data
  set.seed(123)
  n <- 50
  x <- seq(0, 1, length.out = n)
  y <- 2 + rnorm(n, 0, 0.1)
  data <- data.frame(x = x, y = y)
  
  # Fit GAM
  gam_model <- mgcv::gam(y ~ s(x), data = data)
  
  # Test F-test
  result <- test_gam_f_test(gam_model)
  
  # Check structure
  expect_type(result, "list")
  expect_true("method" %in% names(result))
  expect_true("f_statistics" %in% names(result))
  expect_true("p_values" %in% names(result))
  expect_true("is_flat" %in% names(result))
})

test_that("interpret_flatness_results works correctly", {
  # Create mock results
  mock_results <- list(
    test1 = list(is_flat = TRUE),
    test2 = list(is_flat = TRUE),
    test3 = list(is_flat = FALSE)
  )
  
  interpretation <- interpret_flatness_results(mock_results)
  expect_type(interpretation, "character")
  expect_true(nchar(interpretation) > 0)
  
  # Test with all flat
  all_flat <- list(
    test1 = list(is_flat = TRUE),
    test2 = list(is_flat = TRUE)
  )
  
  interpretation_flat <- interpret_flatness_results(all_flat)
  expect_true(grepl("Strong evidence for flatness", interpretation_flat))
  
  # Test with none flat
  none_flat <- list(
    test1 = list(is_flat = FALSE),
    test2 = list(is_flat = FALSE)
  )
  
  interpretation_nonflat <- interpret_flatness_results(none_flat)
  expect_true(grepl("Strong evidence against flatness", interpretation_nonflat))
})
