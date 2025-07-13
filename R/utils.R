#' Transform Data from Single-Cell Experiments
#'
#' Transform data from single-cell experiment objects into a format suitable for analysis.
#' This is a placeholder function that needs to be implemented based on your specific data structure.
#'
#' @param sce1 Single-cell experiment object 1 (reference)
#' @param sce2 Single-cell experiment object 2 (target)
#' @param gene Character string specifying the gene name
#' @param assay Character string specifying the assay to use
#' @param timeCol Character string specifying the column name for time
#' @return Data frame with time, reference, and target columns
#' @export
transform_data <- function(sce1, sce2, gene, peak = NULL,timeCol = "pseudotime",assay = "log_counts", assay2 = NULL) {
  pseudotime = sce1[[timeCol]]
  sce1Gene = sce1[gene]
  if (is.null(assay2)) {
    assay2 = assay
  }
  gene1 = as.numeric(sce1Gene@assays@data[[assay]])
  if (is.null(peak)) {
    region = gene
  }
  else {
    region = peak
  }
  sce2Gene = sce2[region]
  gene2 = as.numeric(sce2Gene@assays@data[[assay2]])
  df = data.frame(time = pseudotime, reference = gene1, target = gene2)
  return(df)
} 