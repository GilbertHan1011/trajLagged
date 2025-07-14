#' GAM Flatness Testing Examples
#' 
#' This file demonstrates how to use the GAM flatness testing methods
#' with both simulated data and real trajectory data.

library(mgcv)
library(trajLagged)

# ============================================================================
# Example 1: Testing with Simulated Data
# ============================================================================

#' Create simulated data for testing
create_test_data <- function(n = 100, noise_level = 0.1, flat = FALSE) {
  set.seed(123)
  x <- seq(0, 1, length.out = n)
  
  if (flat) {
    # Flat relationship (just noise around a constant)
    y <- 2 + rnorm(n, 0, noise_level)
  } else {
    # Non-linear relationship
    y <- 2 + sin(2 * pi * x) + rnorm(n, 0, noise_level)
  }
  
  return(data.frame(x = x, y = y))
}

# Example 1a: Test a flat GAM
cat("=== Example 1a: Testing a Flat GAM ===\n")
flat_data <- create_test_data(n = 100, flat = TRUE)
flat_gam <- gam(y ~ s(x), data = flat_data)

# Run all flatness tests
flat_results <- test_gam_flatness(flat_gam, method = "all")

# Print results
print("Flat GAM Test Results:")
print(flat_results$overall_interpretation)

# Detailed results for each test
for (test_name in names(flat_results)) {
  if (test_name != "overall_interpretation") {
    cat("\n", test_name, ":\n")
    cat("  P-value:", flat_results[[test_name]]$p_value, "\n")
    cat("  Is flat:", flat_results[[test_name]]$is_flat, "\n")
  }
}

# Example 1b: Test a non-flat GAM
cat("\n=== Example 1b: Testing a Non-Flat GAM ===\n")
nonflat_data <- create_test_data(n = 100, flat = FALSE)
nonflat_gam <- gam(y ~ s(nonflat_data$x), data = nonflat_data)

nonflat_results <- test_gam_flatness(nonflat_gam, method = "all")
print("Non-Flat GAM Test Results:")
print(nonflat_results$overall_interpretation)

# ============================================================================
# Example 2: Integration with trajLagged Pipeline
# ============================================================================

#' Example using the trajLagged pipeline (requires real data)
#' This example shows how to integrate flatness testing with your existing workflow

# Uncomment and modify the following code when you have real data:

# cat("\n=== Example 2: Testing Trajectory GAMs ===\n")
# 
# # Test flatness for a specific gene in your trajectory data
# trajectory_results <- test_trajectory_flatness(
#   sce1 = your_sce1_object,
#   sce2 = your_sce2_object, 
#   gene = "YOUR_GENE_NAME",
#   test_both = TRUE,  # Test both reference and target
#   method = "all",
#   alpha = 0.05
# )
# 
# # Print results for reference trajectory
# cat("Reference trajectory flatness:\n")
# print(trajectory_results$reference$flatness_tests$overall_interpretation)
# 
# # Print results for target trajectory  
# cat("Target trajectory flatness:\n")
# print(trajectory_results$target$flatness_tests$overall_interpretation)
# 
# # Create diagnostic plots
# plot_gam_flatness(
#   trajectory_results$reference$gam_model,
#   trajectory_results$reference$flatness_tests,
#   main_title = "Reference Trajectory GAM Diagnostics"
# )

# ============================================================================
# Example 3: Batch Testing Multiple Genes
# ============================================================================

#' Function to test flatness for multiple genes
test_multiple_genes_flatness <- function(sce1, sce2, gene_list, method = "summary", alpha = 0.05) {
  
  results <- list()
  
  for (gene in gene_list) {
    cat("Testing gene:", gene, "\n")
    
    tryCatch({
      # Test flatness for this gene
      gene_results <- test_trajectory_flatness(
        sce1 = sce1, sce2 = sce2, gene = gene,
        test_both = TRUE, method = method, alpha = alpha
      )
      
      # Extract key results
      results[[gene]] <- list(
        reference_flat = gene_results$reference$flatness_tests[[paste0(method, "_test")]]$is_flat,
        target_flat = gene_results$target$flatness_tests[[paste0(method, "_test")]]$is_flat,
        reference_pval = gene_results$reference$flatness_tests[[paste0(method, "_test")]]$p_value,
        target_pval = gene_results$target$flatness_tests[[paste0(method, "_test")]]$p_value
      )
      
    }, error = function(e) {
      warning("Error testing gene ", gene, ": ", e$message)
      results[[gene]] <- list(
        reference_flat = NA, target_flat = NA,
        reference_pval = NA, target_pval = NA
      )
    })
  }
  
  return(results)
}

# Example usage (uncomment when you have real data):
# gene_list <- c("GENE1", "GENE2", "GENE3")  # Replace with your genes
# batch_results <- test_multiple_genes_flatness(sce1, sce2, gene_list)

# ============================================================================
# Example 4: Method Comparison
# ============================================================================

#' Compare different testing methods on the same data
compare_flatness_methods <- function(gam_model) {
  
  methods <- c("summary", "wald_contrasts", "likelihood_ratio", "f_test")
  comparison_results <- data.frame(
    method = methods,
    p_value = NA,
    is_flat = NA,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(methods)) {
    method <- methods[i]
    result <- test_gam_flatness(gam_model, method = method)
    
    # Extract p-value (handling different result structures)
    if (method %in% names(result)) {
      test_result <- result[[method]]
      if ("p_value" %in% names(test_result)) {
        comparison_results$p_value[i] <- test_result$p_value
        comparison_results$is_flat[i] <- test_result$is_flat
      } else if ("p_values" %in% names(test_result)) {
        # For methods that return multiple p-values, take the minimum
        p_vals <- test_result$p_values
        comparison_results$p_value[i] <- min(p_vals, na.rm = TRUE)
        comparison_results$is_flat[i] <- test_result$is_flat
      }
    }
  }
  
  return(comparison_results)
}

# Example: Compare methods on flat data
cat("\n=== Example 4: Method Comparison ===\n")
flat_comparison <- compare_flatness_methods(flat_gam)
print("Flat GAM - Method Comparison:")
print(flat_comparison)

# Example: Compare methods on non-flat data
nonflat_comparison <- compare_flatness_methods(nonflat_gam)
print("\nNon-Flat GAM - Method Comparison:")
print(nonflat_comparison)

# ============================================================================
# Example 5: Visualization
# ============================================================================

cat("\n=== Example 5: Creating Diagnostic Plots ===\n")

# Create plots for flat GAM
plot_gam_flatness(flat_gam, flat_results, "Flat GAM Diagnostics")

# Create plots for non-flat GAM  
plot_gam_flatness(nonflat_gam, nonflat_results, "Non-Flat GAM Diagnostics")

# ============================================================================
# Example 6: Custom Significance Levels
# ============================================================================

cat("\n=== Example 6: Testing with Different Alpha Levels ===\n")

alpha_levels <- c(0.01, 0.05, 0.10)
for (alpha in alpha_levels) {
  result <- test_gam_flatness(flat_gam, method = "summary", alpha = alpha)
  cat("Alpha =", alpha, "- Is flat:", result$summary_test$is_flat, "\n")
}

# ============================================================================
# Utility Functions for Results Analysis
# ============================================================================

#' Summarize flatness test results across multiple genes
summarize_flatness_results <- function(batch_results) {
  
  # Convert to data frame for easier analysis
  df <- do.call(rbind, lapply(names(batch_results), function(gene) {
    result <- batch_results[[gene]]
    data.frame(
      gene = gene,
      ref_flat = result$reference_flat,
      target_flat = result$target_flat,
      ref_pval = result$reference_pval,
      target_pval = result$target_pval,
      stringsAsFactors = FALSE
    )
  }))
  
  # Summary statistics
  summary_stats <- list(
    total_genes = nrow(df),
    ref_flat_count = sum(df$ref_flat, na.rm = TRUE),
    target_flat_count = sum(df$target_flat, na.rm = TRUE),
    ref_flat_prop = mean(df$ref_flat, na.rm = TRUE),
    target_flat_prop = mean(df$target_flat, na.rm = TRUE),
    both_flat_count = sum(df$ref_flat & df$target_flat, na.rm = TRUE),
    neither_flat_count = sum(!df$ref_flat & !df$target_flat, na.rm = TRUE)
  )
  
  return(list(
    detailed_results = df,
    summary = summary_stats
  ))
}

cat("\n=== GAM Flatness Testing Examples Complete ===\n")
cat("Use these examples as templates for your own analysis.\n")
cat("Remember to replace simulated data with your actual trajectory data.\n")
