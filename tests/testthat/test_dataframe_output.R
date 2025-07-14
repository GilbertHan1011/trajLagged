#!/usr/bin/env Rscript

#' Test the new dataframe output functionality for GAM flatness testing

# Load required libraries
library(mgcv)

# Source the flatness testing functions
source("R/gam_flatness_tests.R")

cat("=== Testing Dataframe Output for GAM Flatness Testing ===\n\n")

# ============================================================================
# Example 1: Flat GAM - Dataframe Output
# ============================================================================

cat("1. Testing FLAT GAM with dataframe output\n")
cat("------------------------------------------\n")

set.seed(123)
n <- 100
x_flat <- seq(0, 1, length.out = n)
y_flat <- 2 + rnorm(n, 0, 0.1)  # Constant with small noise
data_flat <- data.frame(x = x_flat, y = y_flat)

# Fit GAM
gam_flat <- gam(y ~ s(x), data = data_flat)

# Test flatness with dataframe output
result_flat_df <- test_gam_flatness(gam_flat, method = "all", output_format = "dataframe")

cat("Dataframe structure:\n")
print(str(result_flat_df))

cat("\nDataframe content:\n")
print(result_flat_df)

# ============================================================================
# Example 2: Non-flat GAM - Dataframe Output
# ============================================================================

cat("\n\n2. Testing NON-FLAT GAM with dataframe output\n")
cat("----------------------------------------------\n")

set.seed(123)
x_nonflat <- seq(0, 2*pi, length.out = n)
y_nonflat <- sin(x_nonflat) + rnorm(n, 0, 0.1)  # Clear sinusoidal pattern
data_nonflat <- data.frame(x = x_nonflat, y = y_nonflat)

# Fit GAM
gam_nonflat <- gam(y ~ s(x), data = data_nonflat)

# Test flatness with dataframe output
result_nonflat_df <- test_gam_flatness(gam_nonflat, method = "all", output_format = "dataframe")

cat("Non-flat GAM results:\n")
print(result_nonflat_df)

# ============================================================================
# Example 3: Compare List vs Dataframe Output
# ============================================================================

cat("\n\n3. Comparing List vs Dataframe Output\n")
cat("-------------------------------------\n")

# List output
result_list <- test_gam_flatness(gam_nonflat, method = "all", output_format = "list")
result_df <- test_gam_flatness(gam_nonflat, method = "all", output_format = "dataframe")

cat("List output structure:\n")
print(names(result_list))

cat("\nDataframe columns:\n")
print(names(result_df))

cat("\nKey statistical values comparison:\n")
cat("Summary p-value (list):", min(result_list$summary_test$p_values), "\n")
cat("Summary p-value (df)  :", result_df$summary_p_value, "\n")

cat("Wald statistic (list) :", result_list$wald_contrasts$wald_statistic, "\n")
cat("Wald statistic (df)   :", result_df$wald_statistic, "\n")

cat("LR statistic (list)   :", result_list$likelihood_ratio$lr_statistic, "\n")
cat("LR statistic (df)     :", result_df$lr_statistic, "\n")

# ============================================================================
# Example 4: Batch Processing with Dataframe Output
# ============================================================================

cat("\n\n4. Batch Processing Example\n")
cat("---------------------------\n")

# Create multiple test datasets
datasets <- list(
  flat1 = data.frame(x = seq(0, 1, length.out = 50), 
                     y = 2 + rnorm(50, 0, 0.1)),
  flat2 = data.frame(x = seq(0, 1, length.out = 50), 
                     y = 3 + rnorm(50, 0, 0.15)),
  nonflat1 = data.frame(x = seq(0, 2*pi, length.out = 50), 
                        y = sin(seq(0, 2*pi, length.out = 50)) + rnorm(50, 0, 0.1)),
  nonflat2 = data.frame(x = seq(0, 1, length.out = 50), 
                        y = seq(0, 1, length.out = 50)^2 + rnorm(50, 0, 0.1))
)

# Process all datasets and combine results
batch_results <- do.call(rbind, lapply(names(datasets), function(name) {
  data <- datasets[[name]]
  gam_model <- gam(y ~ s(x), data = data)
  
  result_df <- test_gam_flatness(gam_model, method = "all", output_format = "dataframe")
  result_df$dataset_name <- name
  
  return(result_df)
}))

cat("Batch processing results:\n")
print(batch_results[, c("dataset_name", "summary_p_value", "wald_p_value", 
                        "lr_p_value", "f_p_value", "min_p_value", 
                        "n_tests_flat", "overall_interpretation")])

# ============================================================================
# Example 5: Summary Statistics
# ============================================================================

cat("\n\n5. Summary Statistics from Dataframe\n")
cat("------------------------------------\n")

# Calculate summary statistics
cat("Summary of minimum p-values:\n")
print(summary(batch_results$min_p_value))

cat("\nProportion of tests suggesting flatness by dataset:\n")
print(batch_results[, c("dataset_name", "prop_tests_flat")])

cat("\nDatasets classified as flat (all tests agree):\n")
flat_datasets <- batch_results[batch_results$n_tests_flat == batch_results$n_tests_total, "dataset_name"]
print(flat_datasets)

cat("\nDatasets classified as non-flat (all tests agree):\n")
nonflat_datasets <- batch_results[batch_results$n_tests_flat == 0, "dataset_name"]
print(nonflat_datasets)

# ============================================================================
# Example 6: Easy Filtering and Analysis
# ============================================================================

cat("\n\n6. Easy Filtering with Dataframe Format\n")
cat("---------------------------------------\n")

# Filter for significant results (any test with p < 0.05)
significant_results <- batch_results[batch_results$min_p_value < 0.05, ]
cat("Number of significant results:", nrow(significant_results), "\n")

# Filter for strong evidence against flatness
strong_nonflat <- batch_results[batch_results$n_tests_flat == 0, ]
cat("Number with strong evidence against flatness:", nrow(strong_nonflat), "\n")

# Show the most significant result
most_significant <- batch_results[which.min(batch_results$min_p_value), ]
cat("Most significant result:\n")
cat("  Dataset:", most_significant$dataset_name, "\n")
cat("  Min p-value:", most_significant$min_p_value, "\n")
cat("  Interpretation:", most_significant$overall_interpretation, "\n")

cat("\n=== Dataframe Output Testing Complete ===\n")
cat("The dataframe format provides:\n")
cat("- Easy batch processing and combination of results\n")
cat("- Simple filtering and analysis\n")
cat("- All statistical values in a structured format\n")
cat("- Compatibility with standard R data analysis workflows\n")
