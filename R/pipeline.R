#' Complete GAM Pipeline
#'
#' Run the complete GAM preprocessing and analysis pipeline.
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param gene Character string specifying the gene name
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @return Numeric value representing the analysis result
#' @export
pipeline_gam <- function(sce1, sce2, gene, binned = FALSE) {
  # Use the new preprocessing function with GAM
  df_binned_gam <- preprocess_gam(sce1, sce2, gene, binned = binned)
  
  # Use the updated get_value function
  get_value_gam(df_binned_gam)
}

#' Pipeline from Data Frame
#'
#' Run the GAM pipeline starting from a data frame.
#'
#' @param df Data frame containing time, reference, and target columns
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @return Numeric value representing the analysis result
#' @export
pipeline_from_df <- function(df, binned = FALSE) {
  df_binned_gam <- preprocess_gam_from_df(df, binned = binned)
  get_value_gam(df_binned_gam)
} 