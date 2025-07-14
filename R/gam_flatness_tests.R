#' Statistical Tests for GAM Flatness
#'
#' This module provides comprehensive methods for testing whether a GAM fitted 
#' with the mgcv package represents a flat (non-significant) relationship.
#'
#' @importFrom mgcv gam s
#' @importFrom stats vcov coef predict pchisq qchisq pf qf
#' @importFrom graphics plot

#' Test GAM Flatness Using Multiple Methods
#'
#' Comprehensive function that applies multiple statistical tests to determine
#' if a GAM represents a flat relationship.
#'
#' @param gam_model A fitted GAM model from mgcv package
#' @param method Character vector specifying which tests to perform. Options:
#'   "summary", "wald_contrasts", "likelihood_ratio", "f_test", "all"
#' @param n_points Integer, number of points for contrast testing (default: 50)
#' @param alpha Numeric, significance level (default: 0.05)
#' @param output_format Character, either "list" (default) or "dataframe"
#' @return List containing test results and interpretations, or dataframe if output_format="dataframe"
#' @export
test_gam_flatness <- function(gam_model, method = "all", n_points = 50, alpha = 0.05, output_format = "list") {

  # Validate inputs
  if (!inherits(gam_model, "gam")) {
    stop("gam_model must be a fitted GAM object from mgcv package")
  }

  valid_methods <- c("summary", "wald_contrasts", "likelihood_ratio", "f_test", "all")
  if (!all(method %in% valid_methods)) {
    stop("method must be one of: ", paste(valid_methods, collapse = ", "))
  }

  valid_formats <- c("list", "dataframe")
  if (!output_format %in% valid_formats) {
    stop("output_format must be one of: ", paste(valid_formats, collapse = ", "))
  }

  if ("all" %in% method) {
    method <- c("summary", "wald_contrasts", "likelihood_ratio", "f_test")
  }

  results <- list()

  # Method 1: Summary-based test (approximate F-test)
  if ("summary" %in% method) {
    results$summary_test <- test_gam_summary(gam_model, alpha)
  }

  # Method 2: Advanced Wald Test with Contrasts
  if ("wald_contrasts" %in% method) {
    results$wald_contrasts <- test_gam_wald_contrasts(gam_model, n_points, alpha)
  }

  # Method 3: Likelihood Ratio Test
  if ("likelihood_ratio" %in% method) {
    results$likelihood_ratio <- test_gam_likelihood_ratio(gam_model, alpha)
  }

  # Method 4: F-test for smooth terms
  if ("f_test" %in% method) {
    results$f_test <- test_gam_f_test(gam_model, alpha)
  }

  # Overall interpretation
  results$overall_interpretation <- interpret_flatness_results(results, alpha)

  # Return as dataframe if requested
  if (output_format == "dataframe") {
    return(convert_results_to_dataframe(results, alpha))
  }

  return(results)
}

#' Method 1: Summary-based Test
#'
#' Extract p-values from GAM summary for smooth terms
#'
#' @param gam_model Fitted GAM model
#' @param alpha Significance level
#' @return List with test results
test_gam_summary <- function(gam_model, alpha = 0.05) {
  
  gam_summary <- summary(gam_model)
  
  # Extract smooth terms information
  smooth_terms <- gam_summary$s.table
  
  if (is.null(smooth_terms) || nrow(smooth_terms) == 0) {
    return(list(
      method = "Summary Test",
      p_values = NA,
      is_flat = NA,
      interpretation = "No smooth terms found in model"
    ))
  }
  
  # Get p-values for smooth terms
  p_values <- smooth_terms[, "p-value"]
  names(p_values) <- rownames(smooth_terms)
  
  # Test for flatness (non-significance)
  is_flat <- all(p_values > alpha)
  
  return(list(
    method = "Summary Test (Approximate F-test)",
    p_values = p_values,
    is_flat = is_flat,
    alpha = alpha,
    interpretation = ifelse(is_flat, 
                           "GAM appears flat (all smooth terms non-significant)",
                           "GAM shows significant non-linear relationship")
  ))
}

#' Method 2: Advanced Wald Test with Contrasts
#'
#' Manual Wald test using contrasts between predicted values
#'
#' @param gam_model Fitted GAM model
#' @param n_points Number of points for contrast testing
#' @param alpha Significance level
#' @return List with test results
test_gam_wald_contrasts <- function(gam_model, n_points = 50, alpha = 0.05) {
  
  tryCatch({
    # --- FIX START ---
    # Correctly check for the existence of smooth terms
    if (is.null(gam_model$smooth) || length(gam_model$smooth) == 0) {
      return(list(
        method = "Wald Contrasts Test",
        wald_statistic = NA,
        p_value = NA,
        is_flat = NA,
        interpretation = "No smooth terms found in the model for contrast testing."
      ))
    }
    # --- FIX END ---
    
    # Extract coefficients and covariance matrix
    beta <- coef(gam_model)
    Vp <- vcov(gam_model)
    
    # Create prediction data (assumes the first smooth term is the target)
    predictor_name <- gam_model$smooth[[1]]$term
    original_data <- gam_model$model
    
    # Handle cases where the predictor might not be in the model frame (e.g., transformations)
    if (!predictor_name %in% names(original_data)) {
        stop(paste("Predictor variable '", predictor_name, "' not found in model data frame.", sep=""))
    }
    predictor_range <- range(original_data[[predictor_name]], na.rm = TRUE)
    
    new_data <- data.frame(
      x = seq(predictor_range[1], predictor_range[2], length.out = n_points)
    )
    names(new_data) <- predictor_name
    
    # Get linear predictor matrix
    lp_matrix <- predict(gam_model, newdata = new_data, type = "lpmatrix")
    
    # Create contrast matrix (consecutive differences)
    n_comparisons <- nrow(lp_matrix) - 1
    L <- matrix(0, nrow = n_comparisons, ncol = ncol(lp_matrix))
    for (i in 1:n_comparisons) {
      L[i, ] <- lp_matrix[i + 1, ] - lp_matrix[i, ]
    }
    
    # Calculate Wald statistic
    L_beta <- L %*% beta
    L_Vp_Lt <- L %*% Vp %*% t(L)
    
    # Use generalized inverse for stability instead of checking determinant
    L_Vp_Lt_inv <- MASS::ginv(L_Vp_Lt)
    wald_stat <- as.numeric(t(L_beta) %*% L_Vp_Lt_inv %*% L_beta)
    
    # Degrees of freedom
    df <- qr(L)$rank
    
    # Calculate p-value
    p_value <- pchisq(wald_stat, df = df, lower.tail = FALSE)
    
    # Test for flatness
    is_flat <- p_value > alpha
    
    return(list(
      method = "Wald Contrasts Test",
      wald_statistic = wald_stat,
      degrees_freedom = df,
      p_value = p_value,
      is_flat = is_flat,
      alpha = alpha,
      n_contrasts = n_comparisons,
      interpretation = ifelse(is_flat,
                              "Fail to reject null: GAM appears flat.",
                              "Reject null: GAM shows significant variation.")
    ))
    
  }, error = function(e) {
    return(list(
      method = "Wald Contrasts Test",
      wald_statistic = NA,
      p_value = NA,
      is_flat = NA,
      interpretation = paste("Error in Wald test:", e$message)
    ))
  })
}
#' Method 3: Likelihood Ratio Test
#'
#' Compare GAM with null model (intercept only)
#'
#' @param gam_model Fitted GAM model
#' @param alpha Significance level
#' @return List with test results
test_gam_likelihood_ratio <- function(gam_model, alpha = 0.05) {
  
  tryCatch({
    # Fit null model (intercept only)
    response_var <- as.character(gam_model$formula)[2]
    null_formula <- as.formula(paste(response_var, "~ 1"))
    null_model <- gam(null_formula, data = gam_model$model)
    
    # Calculate likelihood ratio statistic
    lr_stat <- 2 * (logLik(gam_model) - logLik(null_model))
    
    # Degrees of freedom difference
    df_diff <- gam_model$df.residual - null_model$df.residual
    df_diff <- abs(df_diff)  # Ensure positive
    
    # Calculate p-value
    p_value <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
    
    # Test for flatness
    is_flat <- p_value > alpha
    
    return(list(
      method = "Likelihood Ratio Test",
      lr_statistic = as.numeric(lr_stat),
      degrees_freedom = df_diff,
      p_value = p_value,
      is_flat = is_flat,
      alpha = alpha,
      gam_loglik = as.numeric(logLik(gam_model)),
      null_loglik = as.numeric(logLik(null_model)),
      interpretation = ifelse(is_flat,
                             "GAM not significantly better than null model (appears flat)",
                             "GAM significantly better than null model (not flat)")
    ))
    
  }, error = function(e) {
    return(list(
      method = "Likelihood Ratio Test",
      lr_statistic = NA,
      p_value = NA,
      is_flat = NA,
      interpretation = paste("Error in likelihood ratio test:", e$message)
    ))
  })
}

#' Method 4: F-test for Smooth Terms
#'
#' Extract F-statistics for smooth terms from GAM
#'
#' @param gam_model Fitted GAM model
#' @param alpha Significance level
#' @return List with test results
test_gam_f_test <- function(gam_model, alpha = 0.05) {
  
  gam_summary <- summary(gam_model)
  smooth_terms <- gam_summary$s.table
  
  if (is.null(smooth_terms) || nrow(smooth_terms) == 0) {
    return(list(
      method = "F-test for Smooth Terms",
      f_statistics = NA,
      p_values = NA,
      is_flat = NA,
      interpretation = "No smooth terms found in model"
    ))
  }
  
  # Extract F-statistics and p-values
  f_stats <- smooth_terms[, "F"]
  p_values <- smooth_terms[, "p-value"]
  names(f_stats) <- rownames(smooth_terms)
  names(p_values) <- rownames(smooth_terms)
  
  # Test for flatness
  is_flat <- all(p_values > alpha)
  
  return(list(
    method = "F-test for Smooth Terms",
    f_statistics = f_stats,
    p_values = p_values,
    is_flat = is_flat,
    alpha = alpha,
    interpretation = ifelse(is_flat,
                           "All smooth terms non-significant (GAM appears flat)",
                           "Some smooth terms significant (GAM not flat)")
  ))
}

#' Interpret Overall Flatness Results
#'
#' Provide overall interpretation based on multiple test results
#'
#' @param results List of test results
#' @param alpha Significance level
#' @return Character string with interpretation
interpret_flatness_results <- function(results, alpha = 0.05) {

  # Count how many tests suggest flatness
  flatness_votes <- sapply(results, function(x) {
    if (is.list(x) && "is_flat" %in% names(x)) {
      return(x$is_flat)
    }
    return(NA)
  })

  # Remove NA values
  flatness_votes <- flatness_votes[!is.na(flatness_votes)]

  if (length(flatness_votes) == 0) {
    return("No valid test results available for interpretation")
  }

  n_flat <- sum(flatness_votes)
  n_total <- length(flatness_votes)

  if (n_flat == n_total) {
    return(paste("Strong evidence for flatness: All", n_total, "tests suggest the GAM is flat"))
  } else if (n_flat > n_total / 2) {
    return(paste("Moderate evidence for flatness:", n_flat, "of", n_total, "tests suggest flatness"))
  } else if (n_flat == 0) {
    return(paste("Strong evidence against flatness: All", n_total, "tests suggest the GAM is not flat"))
  } else {
    return(paste("Mixed evidence:", n_flat, "of", n_total, "tests suggest flatness"))
  }
}

#' Convert Test Results to Dataframe
#'
#' Convert the list of test results to a single-row dataframe with all statistical values
#'
#' @param results List of test results from test_gam_flatness
#' @param alpha Significance level used in testing
#' @return Single-row dataframe with all statistical values
convert_results_to_dataframe <- function(results, alpha = 0.05) {

  # Initialize dataframe with basic info
  df <- data.frame(
    alpha = alpha,
    stringsAsFactors = FALSE
  )

  # Extract summary test results
  if ("summary_test" %in% names(results)) {
    summary_res <- results$summary_test
    df$summary_p_value <- ifelse(length(summary_res$p_values) > 0,
                                min(summary_res$p_values, na.rm = TRUE), NA)
    df$summary_is_flat <- summary_res$is_flat
  } else {
    df$summary_p_value <- NA
    df$summary_is_flat <- NA
  }

  # Extract Wald contrasts results
  if ("wald_contrasts" %in% names(results)) {
    wald_res <- results$wald_contrasts
    df$wald_statistic <- wald_res$wald_statistic
    df$wald_degrees_freedom <- wald_res$degrees_freedom
    df$wald_p_value <- wald_res$p_value
    df$wald_is_flat <- wald_res$is_flat
    df$wald_n_contrasts <- wald_res$n_contrasts
  } else {
    df$wald_statistic <- NA
    df$wald_degrees_freedom <- NA
    df$wald_p_value <- NA
    df$wald_is_flat <- NA
    df$wald_n_contrasts <- NA
  }

  # Extract likelihood ratio results
  if ("likelihood_ratio" %in% names(results)) {
    lr_res <- results$likelihood_ratio
    df$lr_statistic <- lr_res$lr_statistic
    df$lr_degrees_freedom <- lr_res$degrees_freedom
    df$lr_p_value <- as.numeric(lr_res$p_value)
    df$lr_is_flat <- lr_res$is_flat
    df$lr_gam_loglik <- lr_res$gam_loglik
    df$lr_null_loglik <- lr_res$null_loglik
  } else {
    df$lr_statistic <- NA
    df$lr_degrees_freedom <- NA
    df$lr_p_value <- NA
    df$lr_is_flat <- NA
    df$lr_gam_loglik <- NA
    df$lr_null_loglik <- NA
  }

  # Extract F-test results
  if ("f_test" %in% names(results)) {
    f_res <- results$f_test
    df$f_statistic <- ifelse(length(f_res$f_statistics) > 0,
                            max(f_res$f_statistics, na.rm = TRUE), NA)
    df$f_p_value <- ifelse(length(f_res$p_values) > 0,
                          min(f_res$p_values, na.rm = TRUE), NA)
    df$f_is_flat <- f_res$is_flat
  } else {
    df$f_statistic <- NA
    df$f_p_value <- NA
    df$f_is_flat <- NA
  }

  # Overall results
  if ("overall_interpretation" %in% names(results)) {
    df$overall_interpretation <- results$overall_interpretation

    # Count how many tests suggest flatness
    flatness_cols <- c("summary_is_flat", "wald_is_flat", "lr_is_flat", "f_is_flat")
    available_tests <- flatness_cols[flatness_cols %in% names(df)]

    if (length(available_tests) > 0) {
      flat_votes <- sapply(available_tests, function(col) df[[col]])
      flat_votes <- flat_votes[!is.na(flat_votes)]

      df$n_tests_total <- length(flat_votes)
      df$n_tests_flat <- sum(flat_votes)
      df$prop_tests_flat <- ifelse(length(flat_votes) > 0,
                                  sum(flat_votes) / length(flat_votes), NA)
    } else {
      df$n_tests_total <- 0
      df$n_tests_flat <- 0
      df$prop_tests_flat <- NA
    }
  } else {
    df$overall_interpretation <- "No interpretation available"
    df$n_tests_total <- 0
    df$n_tests_flat <- 0
    df$prop_tests_flat <- NA
  }

  # Add minimum p-value across all tests
  p_value_cols <- c("summary_p_value", "wald_p_value", "lr_p_value", "f_p_value")
  available_p_values <- p_value_cols[p_value_cols %in% names(df)]

  if (length(available_p_values) > 0) {
    all_p_values <- sapply(available_p_values, function(col) df[[col]])
    all_p_values <- all_p_values[!is.na(all_p_values)]
    df$min_p_value <- ifelse(length(all_p_values) > 0, min(all_p_values), NA)
  } else {
    df$min_p_value <- NA
  }

  return(df)
}

#' Test Flatness for GAMs in trajLagged Pipeline
#'
#' Convenience function to test flatness of GAMs fitted in the trajLagged preprocessing pipeline
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param gene Character string specifying the gene name
#' @param peak Optional peak information
#' @param assay Assay name for sce1 (default: "log_counts")
#' @param assay2 Assay name for sce2 (default: NULL, uses same as assay)
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @param test_both Logical, whether to test both reference and target GAMs (default: TRUE)
#' @param method Character vector specifying which tests to perform (default: "all")
#' @param n_points Integer, number of points for contrast testing (default: 50)
#' @param alpha Numeric, significance level (default: 0.05)
#' @param output_format Character, either "list" (default) or "dataframe"
#' @return List containing flatness test results for reference and/or target GAMs, or dataframe if output_format="dataframe"
#' @export
test_trajectory_flatness <- function(sce1, sce2, gene, peak = NULL, assay = "log_counts", assay2 = NULL,
                                   binned = FALSE, test_both = TRUE, method = "all",
                                   n_points = 50, alpha = 0.05, output_format = "list") {

  # Get the data using existing pipeline
  df <- transform_data(sce1, sce2, gene, peak, assay = assay, assay2 = assay2)

  if (binned) {
    df_binned <- df %>%
      mutate(time_bin = cut(time, breaks = 100, labels = FALSE)) %>%
      group_by(time_bin) %>%
      summarize(
        mean_time = mean(time),
        mean_reference = mean(reference),
        mean_target = mean(target)
      ) %>%
      ungroup()
  } else {
    df_binned <- df
    colnames(df_binned) <- c("mean_time", "mean_reference", "mean_target")
  }

  # Fit GAMs
  gam_ref <- gam(mean_reference ~ s(mean_time), data = df_binned)
  gam_tar <- gam(mean_target ~ s(mean_time), data = df_binned)

  results <- list()

  # Test reference GAM
  ref_flatness <- test_gam_flatness(gam_ref, method = method, n_points = n_points,
                                   alpha = alpha, output_format = output_format)

  results$reference <- list(
    gene = gene,
    model_type = "reference",
    gam_model = gam_ref,
    flatness_tests = ref_flatness
  )

  # Test target GAM if requested
  if (test_both) {
    tar_flatness <- test_gam_flatness(gam_tar, method = method, n_points = n_points,
                                     alpha = alpha, output_format = output_format)

    results$target <- list(
      gene = gene,
      model_type = "target",
      gam_model = gam_tar,
      flatness_tests = tar_flatness
    )
  }

  # If dataframe output requested, combine results
  if (output_format == "dataframe") {
    ref_df <- ref_flatness
    ref_df$gene <- gene
    ref_df$model_type <- "reference"

    if (test_both) {
      tar_df <- tar_flatness
      tar_df$gene <- gene
      tar_df$model_type <- "target"

      # Combine reference and target results
      combined_df <- rbind(ref_df, tar_df)
      return(combined_df)
    } else {
      return(ref_df)
    }
  }

  return(results)
}

#' Visualize GAM Flatness Test Results
#'
#' Create diagnostic plots for GAM flatness testing
#'
#' @param gam_model Fitted GAM model
#' @param flatness_results Results from test_gam_flatness()
#' @param main_title Main title for the plot
#' @return Invisible NULL (creates plots)
#' @importFrom graphics par layout plot abline text mtext
#' @export
plot_gam_flatness <- function(gam_model, flatness_results = NULL, main_title = "GAM Flatness Diagnostic") {

  # Set up plotting layout
  old_par <- par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  on.exit(par(old_par))

  # Plot 1: GAM smooth with confidence intervals
  plot(gam_model, main = "GAM Smooth Term",
       shade = TRUE, residuals = TRUE, pch = 1, cex = 0.8)
  abline(h = 0, col = "red", lty = 2)

  # Plot 2: Residuals vs fitted
  fitted_vals <- fitted(gam_model)
  residuals_vals <- residuals(gam_model)
  plot(fitted_vals, residuals_vals,
       xlab = "Fitted Values", ylab = "Residuals",
       main = "Residuals vs Fitted", pch = 16, cex = 0.8)
  abline(h = 0, col = "red", lty = 2)

  # Plot 3: QQ plot of residuals
  qqnorm(residuals_vals, main = "Normal Q-Q Plot", pch = 16, cex = 0.8)
  qqline(residuals_vals, col = "red")

  # Plot 4: Test results summary
  if (!is.null(flatness_results)) {
    plot.new()

    # Extract p-values from different tests
    test_names <- c()
    p_values <- c()

    if ("summary_test" %in% names(flatness_results)) {
      test_names <- c(test_names, "Summary")
      p_val <- flatness_results$summary_test$p_values
      p_values <- c(p_values, ifelse(length(p_val) > 0, min(p_val), NA))
    }

    if ("wald_contrasts" %in% names(flatness_results)) {
      test_names <- c(test_names, "Wald")
      p_values <- c(p_values, flatness_results$wald_contrasts$p_value)
    }

    if ("likelihood_ratio" %in% names(flatness_results)) {
      test_names <- c(test_names, "LR")
      p_values <- c(p_values, flatness_results$likelihood_ratio$p_value)
    }

    if ("f_test" %in% names(flatness_results)) {
      test_names <- c(test_names, "F-test")
      p_val <- flatness_results$f_test$p_values
      p_values <- c(p_values, ifelse(length(p_val) > 0, min(p_val), NA))
    }

    # Create text summary
    text(0.1, 0.9, "Flatness Test Results:", cex = 1.2, font = 2)

    y_pos <- 0.8
    for (i in seq_along(test_names)) {
      if (!is.na(p_values[i])) {
        result_text <- sprintf("%s: p = %.4f", test_names[i], p_values[i])
        color <- ifelse(p_values[i] > 0.05, "blue", "red")
        text(0.1, y_pos, result_text, cex = 1, col = color)
        y_pos <- y_pos - 0.1
      }
    }

    # Overall interpretation
    if ("overall_interpretation" %in% names(flatness_results)) {
      text(0.1, y_pos - 0.1, "Overall:", cex = 1.1, font = 2)
      text(0.1, y_pos - 0.2, flatness_results$overall_interpretation,
           cex = 0.9, col = "darkgreen")
    }
  }

  # Add main title
  mtext(main_title, outer = TRUE, cex = 1.3, font = 2, line = -1.5)

  invisible(NULL)
}
