#' GAM-based Preprocessing for Single-Cell Trajectory Data
#'
#' This function performs GAM-based smoothing and preprocessing of single-cell
#' trajectory data with optional binning and appropriate error distributions.
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param gene Character string specifying the gene name
#' @param peak Optional peak information
#' @param assay Assay name for sce1 (default: "log_counts")
#' @param assay2 Assay name for sce2 (default: NULL, uses same as assay)
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @param family GAM family specification (default: "auto" for automatic selection)
#'   Options: "auto", "gaussian", "nb", "poisson", "Gamma", "tw"
#' @return A data frame with smoothed time series data
#' @importFrom mgcv gam s nb tw
#' @importFrom dplyr %>% mutate group_by summarize ungroup
#' @importFrom stats predict gaussian poisson Gamma
#' @export
preprocess_gam <- function(sce1, sce2, gene, peak = NULL, assay = "log_counts",
                          assay2 = NULL, binned = FALSE, family = "nb") {
  # Original data transformation and binning
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
      ungroup() # Ungroup for cleaner subsequent steps
  } else {
    df_binned <- df
    colnames(df_binned) <- c("mean_time", "mean_reference", "mean_target")
    df_binned$time_bin <- df_binned$mean_time
    df_binned <- df_binned[c("time_bin", "mean_time", "mean_reference", "mean_target")]
  }
  
  #--- GAM SMOOTHING STEP ---#
  # Select appropriate family/link function for gene expression data
  gam_family <- select_gam_family(df_binned, family)

  # Smooth the reference and target data against time with appropriate family
  gam_ref <- gam(mean_reference ~ s(mean_time), data = df_binned, family = gam_family)
  gam_tar <- gam(mean_target ~ s(mean_time), data = df_binned, family = gam_family)
  
  prediction_grid <- data.frame(
    mean_time = seq(min(df_binned$mean_time), max(df_binned$mean_time), length.out = 100)
  )
  
  # Predict the smoothed values onto this new grid
  ref_predicted <- predict(gam_ref, newdata = prediction_grid)
  tar_predicted <- predict(gam_tar, newdata = prediction_grid)
  
  # --- Create the final output data frame from the predicted values ---
  df_smoothed <- data.frame(
    time = prediction_grid$mean_time,
    reference_smoothed = ref_predicted,
    target_smoothed = tar_predicted,
    sum_weight = 1 - prediction_grid$mean_time
  )
  
  return(df_smoothed)
}

#' Select Appropriate GAM Family for Gene Expression Data
#'
#' Automatically selects the most appropriate error distribution and link function
#' for gene expression data based on data characteristics.
#'
#' @param data Data frame containing expression values
#' @param family_choice Character string specifying family choice
#' @return GAM family object
#' @importFrom mgcv nb tw
#' @importFrom stats gaussian poisson Gamma
select_gam_family <- function(data, family_choice = "auto") {

  if (family_choice != "auto") {
    # User specified family
    return(switch(family_choice,
                  "gaussian" = gaussian(),
                  "nb" = nb(),
                  "poisson" = poisson(),
                  "Gamma" = Gamma(link = "log"),
                  "tw" = tw(),
                  gaussian()  # default fallback
    ))
  }

  # Automatic family selection based on data characteristics
  ref_values <- data$mean_reference
  tar_values <- data$mean_target
  all_values <- c(ref_values, tar_values)

  # Remove any NA or infinite values
  all_values <- all_values[is.finite(all_values)]

  if (length(all_values) == 0) {
    warning("No finite values found, using Gaussian family")
    return(gaussian())
  }

  # Check for negative values
  has_negative <- any(all_values < 0)

  # Calculate basic statistics
  mean_val <- mean(all_values)
  var_val <- var(all_values)
  min_val <- min(all_values)

  # Decision logic for family selection
  if (has_negative) {
    # If negative values present, use Gaussian
    return(gaussian())
  } else if (min_val >= 0 && var_val > mean_val * 1.5) {
    # Overdispersed count-like data: use Negative Binomial
    return(nb())
  } else if (min_val >= 0 && all(all_values == round(all_values))) {
    # Count data with variance ≈ mean: use Poisson
    return(poisson())
  } else if (min_val > 0 && var_val > mean_val) {
    # Positive continuous data with overdispersion: use Gamma
    return(Gamma(link = "log"))
  } else {
    # Default to Gaussian for other cases
    return(gaussian())
  }
}

#' GAM-based Preprocessing from Data Frame
#'
#' Process trajectory data directly from a data frame with GAM smoothing and
#' appropriate error distributions.
#'
#' @param df Data frame containing time, reference, and target columns
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @param family GAM family specification (default: "auto" for automatic selection)
#'   Options: "auto", "gaussian", "nb", "poisson", "Gamma", "tw"
#' @return A data frame with smoothed time series data
#' @importFrom mgcv gam s nb tw
#' @importFrom dplyr %>% mutate group_by summarize ungroup
#' @importFrom stats predict gaussian poisson Gamma
#' @export
preprocess_gam_from_df <- function(df, binned = FALSE, family = "nb") {
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
    df_binned$time_bin <- df_binned$mean_time
    df_binned <- df_binned[c("time_bin", "mean_time", "mean_reference", "mean_target")]
  }
  
  #--- GAM SMOOTHING STEP ---#
  # Select appropriate family/link function for gene expression data
  gam_family <- select_gam_family(df_binned, family)

  # Smooth the reference and target data against time with appropriate family
  gam_ref <- gam(mean_reference ~ s(mean_time), data = df_binned, family = gam_family)
  gam_tar <- gam(mean_target ~ s(mean_time), data = df_binned, family = gam_family)
  
  prediction_grid <- data.frame(
    mean_time = seq(min(df_binned$mean_time), max(df_binned$mean_time), length.out = 100)
  )
  
  ref_predicted <- predict(gam_ref, newdata = prediction_grid)
  tar_predicted <- predict(gam_tar, newdata = prediction_grid)
  
  df_smoothed <- data.frame(
    time = prediction_grid$mean_time,
    reference_smoothed = ref_predicted,
    target_smoothed = tar_predicted,
    sum_weight = 1 - prediction_grid$mean_time
  )
  
  return(df_smoothed)
} 