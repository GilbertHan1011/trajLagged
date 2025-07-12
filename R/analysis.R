#' Weighted Sum Calculation
#'
#' Calculate a weighted sum with normalization.
#'
#' @param x Numeric vector of values
#' @param weight Numeric vector of weights
#' @return Numeric value representing the weighted sum
#' @export
sumWeight <- function(x, weight) {
  x <- unlist(x)
  x <- (x - min(x)) / (max(x) - min(x))
  weight <- unlist(weight)
  sum((x / sum(x)) * weight)
}

#' Get Value from GAM Analysis
#'
#' Calculate the difference between target and reference weighted sums.
#'
#' @param df_gam Data frame with GAM smoothed results
#' @return Numeric value representing the difference
#' @export
get_value_gam <- function(df_gam) {
  # Use the new smoothed column names
  ref <- sumWeight(df_gam$reference_smoothed, df_gam$sum_weight)
  target <- sumWeight(df_gam$target_smoothed, df_gam$sum_weight)
  
  return(target - ref)
}

#' Calculate Pseudotime Density Weights
#'
#' Calculate weights based on pseudotime density for bootstrap sampling.
#'
#' @param sce1 Single-cell experiment object containing pseudotime information
#' @return Numeric vector of weights
#' @importFrom stats density approx
#' @export
calculate_pseudotime_weights <- function(sce1) {
  pseudotime <- sce1$pseudotime
  dens <- density(pseudotime, n = 2048) # Use high resolution for accuracy
  interp_dens <- approx(dens$x, dens$y, xout = pseudotime)$y
  weights <- 1 / (interp_dens + 1e-6)
  return(weights)
} 