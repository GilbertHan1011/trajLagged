#' trajLagged: GAM-based Preprocessing and Bootstrap Analysis for Single-Cell Trajectory Data
#'
#' This package provides functions for GAM-based smoothing and preprocessing of single-cell 
#' trajectory data, along with bootstrap analysis methods for statistical inference.
#' 
#' @details
#' The package includes functions for:
#' \itemize{
#'   \item Data transformation from single-cell experiment objects
#'   \item GAM smoothing of trajectory data
#'   \item Weighted sum calculations
#'   \item Bootstrap confidence interval estimation with pseudotime density weighting
#' }
#' 
#' Main functions include:
#' \itemize{
#'   \item \code{\link{preprocess_gam}}: GAM-based preprocessing
#'   \item \code{\link{pipeline_gam}}: Complete analysis pipeline
#'   \item \code{\link{bootstrap_pipeline_gam}}: Bootstrap analysis with statistical inference
#' }
#'
#' @keywords internal
"_PACKAGE" 