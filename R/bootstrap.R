#' Bootstrap Pipeline with GAM Analysis
#'
#' Perform bootstrap analysis with m-out-of-n sampling and pseudotime density weighting.
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param gene Character string specifying the gene name
#' @param peak Optional peak information
#' @param assay Assay name for sce1 (default: "log_counts")
#' @param assay2 Assay name for sce2 (default: NULL, uses same as assay)
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @param n_bootstrap Integer, number of bootstrap iterations (default: 100)
#' @param m_prop Numeric, proportion of data to sample in each bootstrap
#'   (default: 0.8)
#' @param test_flatness Logical, whether to perform GAM flatness testing
#'   (default: FALSE)
#' @param flatness_method Character vector specifying which flatness tests to
#'   perform (default: "summary")
#' @param flatness_alpha Numeric, significance level for flatness testing
#'   (default: 0.05)
#' @param family GAM family specification (default: "auto")
#' @return List containing bootstrap results and statistical measures
#' @importFrom stats quantile sd
#' @export
bootstrap_pipeline_gam <- function(sce1, sce2, gene, peak = NULL,
                                  assay = "log_counts", assay2 = NULL,
                                  binned = FALSE, n_bootstrap = 100,
                                  m_prop = 0.8, test_flatness = FALSE,
                                  flatness_method = "all",
                                  flatness_alpha = 0.05, family = "nb") {
  
  # 1. Run the slow data extraction ONCE before the loop
  df_full <- transform_data(sce1, sce2, gene, peak,assay = assay, assay2 = assay2)
  
  # Calculate pseudotime density weights
  weights <- calculate_pseudotime_weights(sce1)

  # GAM flatness testing on original data (if requested)
  flatness_results <- NULL
  if (test_flatness) {
    flatness_results <- test_trajectory_flatness(
      sce1, sce2, gene, peak, assay, assay2, binned = FALSE,
      test_both = TRUE, method = flatness_method, alpha = flatness_alpha,
      output_format = "dataframe"
    )
  }

  # Calculate m for m-out-of-n bootstrap
  n_total <- nrow(df_full)
  m_size <- round(n_total * m_prop)
  
  # 2. Bootstrap by resampling the DATA FRAME, which is very fast
  bootstrap_results <- replicate(n_bootstrap, {
    # Resample rows from the data frame using pseudotime density weights
    boot_indices <- sample(1:nrow(df_full), size = m_size, replace = TRUE, prob = weights)
    df_boot <- df_full[boot_indices, ]
    
    # Run the FASTER pipeline version that starts from a data frame
    tryCatch({
      pipeline_from_df(df_boot, binned = binned, family = family)
    }, error = function(e) NA)
  })
  
  # Remove failed bootstraps
  bootstrap_results <- bootstrap_results[!is.na(bootstrap_results)]
  
  # Calculate confidence intervals
  ci_lower <- quantile(bootstrap_results, 0.025)
  ci_upper <- quantile(bootstrap_results, 0.975)
  
  # Statistical test: does CI exclude zero?
  excludes_zero <- ci_lower > 0 | ci_upper < 0
  
  result_list <- list(
    bootstrap_values = bootstrap_results,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    excludes_zero = excludes_zero,
    bootstrap_mean = mean(bootstrap_results),
    bootstrap_sd = sd(bootstrap_results),
    p_value_approx = 2 * min(mean(bootstrap_results <= 0),
                            mean(bootstrap_results >= 0)),
    m_size = m_size,
    n_total = n_total
  )

  # Add flatness testing results if performed
  if (test_flatness && !is.null(flatness_results)) {
    result_list$flatness_tests = flatness_results
  }

  return(result_list)
} 