#!/usr/bin/env Rscript

#' Example: GAM Flatness Testing in Bootstrap and Whole Pipeline
#' 
#' This script demonstrates how to use the new GAM flatness testing functionality
#' integrated into the bootstrap and whole gene pipeline functions.

library(mgcv)
library(trajLagged)

cat("=== GAM Flatness Testing in Bootstrap Pipeline ===\n\n")

# ============================================================================
# Example 1: Single Gene Bootstrap with Flatness Testing
# ============================================================================

cat("1. Single Gene Bootstrap Analysis with Flatness Testing\n")
cat("-------------------------------------------------------\n")

# Note: This example uses simulated data. Replace with your actual SCE objects.
# 
# Example usage with real data:
# result_with_flatness <- bootstrap_pipeline_gam(
#   sce1 = your_sce1_object,
#   sce2 = your_sce2_object,
#   gene = "YOUR_GENE_NAME",
#   n_bootstrap = 100,
#   test_flatness = TRUE,
#   flatness_method = "all",
#   flatness_alpha = 0.05
# )

# For demonstration, let's show the structure of what you would get:
cat("Expected structure of bootstrap results with flatness testing:\n")
cat("$bootstrap_values      - Vector of bootstrap results\n")
cat("$ci_lower, $ci_upper   - Confidence intervals\n")
cat("$excludes_zero         - Whether CI excludes zero\n")
cat("$bootstrap_mean        - Mean of bootstrap distribution\n")
cat("$bootstrap_sd          - Standard deviation of bootstrap\n")
cat("$p_value_approx        - Approximate p-value\n")
cat("$m_size, $n_total      - Sample size information\n")
cat("$flatness_tests        - Dataframe with flatness test results\n")

cat("\nFlatness test dataframe contains:\n")
cat("- alpha: Significance level used\n")
cat("- summary_p_value: P-value from summary test\n")
cat("- wald_statistic, wald_p_value: Wald test results\n")
cat("- lr_statistic, lr_p_value: Likelihood ratio test results\n")
cat("- f_statistic, f_p_value: F-test results\n")
cat("- min_p_value: Minimum p-value across all tests\n")
cat("- n_tests_flat: Number of tests suggesting flatness\n")
cat("- overall_interpretation: Text interpretation\n")

# ============================================================================
# Example 2: Whole Gene Pipeline with Flatness Testing
# ============================================================================

cat("\n\n2. Whole Gene Pipeline with Flatness Testing\n")
cat("--------------------------------------------\n")

# Example usage with real data:
# results_with_flatness <- run_whole_gene(
#   sce1 = your_sce1_object,
#   sce2 = your_sce2_object,
#   geneList = c("GENE1", "GENE2", "GENE3"),
#   bootstrap = TRUE,
#   n_bootstrap = 100,
#   test_flatness = TRUE,
#   flatness_method = "summary",  # Use "all" for comprehensive testing
#   flatness_alpha = 0.05,
#   n_cores = 1
# )

cat("Expected structure of whole gene results with flatness testing:\n")
cat("$valueArr              - Named vector of pipeline values\n")
cat("$arrayDf               - Dataframe of bootstrap arrays\n")
cat("$pvalueArr             - Named vector of p-values\n")
cat("$success_rate          - Proportion of successful analyses\n")
cat("$flatness_tests        - Named list of flatness test dataframes\n")

cat("\nAccessing flatness results:\n")
cat("# Get flatness results for a specific gene\n")
cat("gene_flatness <- results_with_flatness$flatness_tests[['GENE1']]\n")
cat("print(gene_flatness)\n")

# ============================================================================
# Example 3: Analyzing Flatness Results
# ============================================================================

cat("\n\n3. Analyzing Flatness Test Results\n")
cat("----------------------------------\n")

cat("Example analysis workflow:\n\n")

cat("# Extract all flatness results into a single dataframe\n")
cat("all_flatness <- do.call(rbind, lapply(names(results_with_flatness$flatness_tests), function(gene) {\n")
cat("  df <- results_with_flatness$flatness_tests[[gene]]\n")
cat("  df$gene <- gene\n")
cat("  return(df)\n")
cat("}))\n\n")

cat("# Identify genes with flat reference trajectories\n")
cat("flat_reference <- all_flatness[all_flatness$model_type == 'reference' & \n")
cat("                              all_flatness$summary_p_value > 0.05, ]\n")
cat("print(flat_reference$gene)\n\n")

cat("# Identify genes with significant non-flat target trajectories\n")
cat("nonflat_target <- all_flatness[all_flatness$model_type == 'target' & \n")
cat("                              all_flatness$min_p_value < 0.05, ]\n")
cat("print(nonflat_target$gene)\n\n")

cat("# Summary statistics\n")
cat("cat('Proportion of flat reference trajectories:', \n")
cat("    mean(all_flatness[all_flatness$model_type == 'reference', 'summary_p_value'] > 0.05))\n")
cat("cat('Proportion of flat target trajectories:', \n")
cat("    mean(all_flatness[all_flatness$model_type == 'target', 'summary_p_value'] > 0.05))\n")

# ============================================================================
# Example 4: Different Testing Strategies
# ============================================================================

cat("\n\n4. Different Testing Strategies\n")
cat("-------------------------------\n")

cat("Strategy 1: Fast screening with summary test only\n")
cat("run_whole_gene(..., test_flatness = TRUE, flatness_method = 'summary')\n\n")

cat("Strategy 2: Comprehensive testing with all methods\n")
cat("run_whole_gene(..., test_flatness = TRUE, flatness_method = 'all')\n\n")

cat("Strategy 3: Custom significance level\n")
cat("run_whole_gene(..., test_flatness = TRUE, flatness_alpha = 0.01)\n\n")

cat("Strategy 4: Non-bootstrap analysis with flatness testing\n")
cat("run_whole_gene(..., bootstrap = FALSE, test_flatness = TRUE)\n")

# ============================================================================
# Example 5: Filtering and Visualization
# ============================================================================

cat("\n\n5. Filtering and Visualization Examples\n")
cat("---------------------------------------\n")

cat("# Filter for genes with strong evidence against flatness\n")
cat("strong_nonflat <- all_flatness[all_flatness$n_tests_flat == 0, ]\n\n")

cat("# Filter for genes with mixed evidence\n")
cat("mixed_evidence <- all_flatness[all_flatness$n_tests_flat > 0 & \n")
cat("                              all_flatness$n_tests_flat < all_flatness$n_tests_total, ]\n\n")

cat("# Create summary plot (requires ggplot2)\n")
cat("# library(ggplot2)\n")
cat("# ggplot(all_flatness, aes(x = model_type, y = -log10(min_p_value))) +\n")
cat("#   geom_boxplot() +\n")
cat("#   geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = 'red') +\n")
cat("#   labs(title = 'GAM Flatness Test Results',\n")
cat("#        y = '-log10(minimum p-value)',\n")
cat("#        x = 'Model Type') +\n")
cat("#   theme_minimal()\n")

# ============================================================================
# Example 6: Integration with Existing Workflow
# ============================================================================

cat("\n\n6. Integration with Existing Workflow\n")
cat("-------------------------------------\n")

cat("# Your existing workflow:\n")
cat("# results <- run_whole_gene(sce1, sce2, geneList, bootstrap = TRUE)\n")
cat("# significant_genes <- names(results$pvalueArr)[results$pvalueArr < 0.05]\n\n")

cat("# Enhanced workflow with flatness testing:\n")
cat("results_enhanced <- run_whole_gene(\n")
cat("  sce1, sce2, geneList, \n")
cat("  bootstrap = TRUE,\n")
cat("  test_flatness = TRUE,\n")
cat("  flatness_method = 'summary'\n")
cat(")\n\n")

cat("# Identify significant AND non-flat genes\n")
cat("significant_genes <- names(results_enhanced$pvalueArr)[results_enhanced$pvalueArr < 0.05]\n")
cat("flatness_df <- do.call(rbind, results_enhanced$flatness_tests)\n")
cat("nonflat_genes <- unique(flatness_df[flatness_df$summary_p_value < 0.05, 'gene'])\n")
cat("final_candidates <- intersect(significant_genes, nonflat_genes)\n")

# ============================================================================
# Performance Considerations
# ============================================================================

cat("\n\n7. Performance Considerations\n")
cat("-----------------------------\n")

cat("- Flatness testing adds computational overhead\n")
cat("- Use flatness_method = 'summary' for fastest results\n")
cat("- Use flatness_method = 'all' for most robust results\n")
cat("- Consider running flatness tests only on significant genes\n")
cat("- Parallel processing (n_cores > 1) recommended for large gene lists\n")

cat("\n=== Example Complete ===\n")
cat("Remember to:\n")
cat("1. Replace simulated data with your actual SCE objects\n")
cat("2. Adjust gene lists to your specific genes of interest\n")
cat("3. Choose appropriate significance levels and testing methods\n")
cat("4. Consider computational resources for large-scale analyses\n")
