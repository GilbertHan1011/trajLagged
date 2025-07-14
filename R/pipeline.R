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
pipeline_gam <- function(sce1, sce2, gene, peak = NULL, binned = FALSE) {
  # Use the new preprocessing function with GAM
  df_binned_gam <- preprocess_gam(sce1, sce2, gene, peak, binned = binned)
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

#' Run the GAM pipeline for a list of genes
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param geneList List of gene names
#' @param peakList List of peak names (optional)
#' @param binned Logical, whether to bin the data (default: FALSE)
#' @param bootstrap Logical, whether to bootstrap the data (default: TRUE)
#' @param n_bootstrap Number of bootstrap samples (default: 100)
#' @param m_prop Proportion of data to sample for bootstrap (default: 0.5)
#' @param n_cores Number of cores to use for parallel processing (default: 1)
#' @param assay Assay name for sce1 (default: "log_counts")
#' @param assay2 Assay name for sce2 (default: NULL, uses same as assay)
#' @param test_flatness Logical, whether to perform GAM flatness testing
#'   (default: FALSE)
#' @param flatness_method Character vector specifying which flatness tests to
#'   perform (default: "all")
#' @param flatness_alpha Numeric, significance level for flatness testing
#'   (default: 0.05)
#' @return List of values, arrays, and p-values for each gene
#' @importFrom parallel mclapply detectCores
#' @export
run_whole_gene <- function(sce1, sce2, geneList, peakList = NULL,
                          binned = FALSE, bootstrap = TRUE, n_bootstrap = 100,
                          m_prop = 0.5, n_cores = 1, assay = "log_counts",
                          assay2 = NULL, test_flatness = FALSE,
                          flatness_method = "all", flatness_alpha = 0.05) {
  
  # Validate inputs
  if (!is.character(geneList) || length(geneList) == 0) {
    stop("geneList must be a non-empty character vector")
  }
  if (!is.null(peakList) && length(peakList) != length(geneList)) {
    stop("peakList must be the same length as geneList")
  }
  # Set up parallel processing
  if (n_cores == 1) {
    n_cores <- 1  # Sequential processing
  } else if (n_cores <= 0) {
    n_cores <- max(1, detectCores() - 1)  # Use all but one core
  }
  
  # Progress tracking function
  progress_msg <- function(i, total) {
    cat(sprintf("Processing gene %d of %d (%s)\n", i, total, geneList[i]))
  }
  
  # Define worker function for processing individual genes
  process_gene <- function(i) {
    gene <- geneList[i]
    if (!is.null(peakList)) {
      peak <- peakList[i]
    }
    else {
      peak <- NULL
    }
    
    # Progress reporting (only in sequential mode to avoid conflicts)
    if (n_cores == 1) {
      progress_msg(i, length(geneList))
    }
    
    tryCatch({
      if (bootstrap) {
        res <- bootstrap_pipeline_gam(sce1, sce2, gene, peak,
                                    assay = assay, assay2 = assay2,
                                    binned = binned,
                                    n_bootstrap = n_bootstrap,
                                    m_prop = m_prop,
                                    test_flatness = test_flatness,
                                    flatness_method = flatness_method,
                                    flatness_alpha = flatness_alpha)
        result <- list(
          value = res$bootstrap_mean,
          bootstrap_values = res$bootstrap_values,
          p_value = res$p_value_approx,
          success = TRUE
        )

        # Add flatness test results if available
        if (test_flatness && !is.null(res$flatness_tests)) {
          result$flatness_tests <- res$flatness_tests
        }

        result
      } else {
        result <- list(
          value = pipeline_gam(sce1, sce2, gene, peak, binned = binned,
                              assay = assay, assay2 = assay2),
          success = TRUE
        )

        # Add flatness test results for non-bootstrap case if requested
        if (test_flatness) {
          flatness_res <- test_trajectory_flatness(
            sce1, sce2, gene, peak, assay, assay2, binned,
            test_both = TRUE, method = flatness_method,
            alpha = flatness_alpha, output_format = "dataframe"
          )
          result$flatness_tests <- flatness_res
        }

        result
      }
    }, error = function(e) {
      warning(sprintf("Error processing gene %s: %s", gene, e$message))
      list(
        value = NA,
        bootstrap_values = if (bootstrap) numeric(0) else NULL,
        p_value = if (bootstrap) NA else NULL,
        flatness_tests = if (test_flatness) NULL else NULL,
        success = FALSE
      )
    })
  }
  
  # Run processing (parallel or sequential)
  if (n_cores > 1) {
    cat(sprintf("Processing %d genes using %d cores...\n", length(geneList), n_cores))
    results <- mclapply(seq_along(geneList), process_gene, mc.cores = n_cores)
  } else {
    results <- lapply(seq_along(geneList), process_gene)
  }
  if (is.null(peakList)) {
    name = geneList
  }
  else {
    name = paste0(geneList, "_", peakList)
  }
  # Extract results
  valueArr <- sapply(results, function(x) x$value)
  names(valueArr) <- name
  
  if (bootstrap) {
    # Extract bootstrap arrays and p-values
    arrayList <- lapply(results, function(x) x$bootstrap_values)
    pvalueArr <- sapply(results, function(x) x$p_value)
    names(pvalueArr) <- name

    # Extract flatness test results if available
    flatness_list <- NULL
    if (test_flatness) {
      flatness_list <- lapply(results, function(x) x$flatness_tests)
      names(flatness_list) <- name
      flatnessDf = do.call(rbind, flatness_list)
    }
    
    # Convert bootstrap results to data frame
    # Handle cases where some genes failed or have different lengths
    valid_arrays <- arrayList[sapply(arrayList, function(x) length(x) > 0)]
    
    if (length(valid_arrays) > 0) {
      max_len <- max(sapply(valid_arrays, length))
      
      # Pad arrays to same length and combine
      arrayMat <- do.call(cbind, lapply(seq_along(geneList), function(i) {
        arr <- arrayList[[i]]
        if (length(arr) == 0) {
          rep(NA, max_len)
        } else {
          c(arr, rep(NA, max_len - length(arr)))
        }
      }))
      
      arrayDf <- as.data.frame(arrayMat)
      colnames(arrayDf) <- name
    } else {
      # All genes failed - create empty data frame
      arrayDf <- data.frame(matrix(NA, nrow = 0, ncol = length(geneList)))
      colnames(arrayDf) <- name
    }
    
    # Report success rate
    success_rate <- mean(sapply(results, function(x) x$success))
    cat(sprintf("Processing completed. Success rate: %.1f%%\n", success_rate * 100))
    
    result_list <- list(
      valueArr = valueArr,
      arrayDf = arrayDf,
      pvalueArr = pvalueArr,
      success_rate = success_rate
    )

    # Add flatness test results if available
    if (test_flatness && !is.null(flatness_list)) {
      result_list$flatness_tests <- flatnessDf
    }

    return(result_list)
  } else {
    success_rate <- mean(sapply(results, function(x) x$success))
    cat(sprintf("Processing completed. Success rate: %.1f%%\n", success_rate * 100))

    result_list <- list(
      valueArr = valueArr,
      success_rate = success_rate
    )

    # Add flatness test results for non-bootstrap case
    if (test_flatness) {
      flatness_list <- lapply(results, function(x) x$flatness_tests)
      names(flatness_list) <- name
      flatnessDf = do.call(rbind, flatness_list)
      result_list$flatness_tests <- flatnessDf
    }

    return(result_list)
  }
}