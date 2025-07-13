#' GAM-based Preprocessing for Single-Cell Trajectory Data
#'
#' This function performs GAM-based smoothing and preprocessing of single-cell 
#' trajectory data with optional binning.
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param gene Character string specifying the gene name
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @return A data frame with smoothed time series data
#' @importFrom mgcv gam s
#' @importFrom dplyr %>% mutate group_by summarize ungroup
#' @importFrom stats predict
#' @export
preprocess_gam <- function(sce1, sce2, gene, peak = NULL, assay = "log_counts", assay2 = NULL, binned = FALSE) {
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
  # Smooth the reference and target data against time
  gam_ref <- gam(mean_reference ~ s(mean_time), data = df_binned)
  gam_tar <- gam(mean_target ~ s(mean_time), data = df_binned)
  
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

#' GAM-based Preprocessing from Data Frame
#'
#' Process trajectory data directly from a data frame with GAM smoothing.
#'
#' @param df Data frame containing time, reference, and target columns
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @return A data frame with smoothed time series data
#' @importFrom mgcv gam s
#' @importFrom dplyr %>% mutate group_by summarize ungroup
#' @importFrom stats predict
#' @export
preprocess_gam_from_df <- function(df, binned = FALSE) {
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
  gam_ref <- gam(mean_reference ~ s(mean_time), data = df_binned)
  gam_tar <- gam(mean_target ~ s(mean_time), data = df_binned)
  
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