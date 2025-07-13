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
transform_data <- function(sce1, sce2, gene,timeCol = "pseudotime",assay = "log_counts") {
  pseudotime = sce1[[timeCol]]
  sce1Gene = sce1[gene]
  sce2Gene = sce2[gene]
  gene1 = sce1[gene]@assays@data[[assay]] %>% as.numeric()
  gene2 = sce2[gene]@assays@data[[assay]] %>% as.numeric()
  df = data.frame(time = pseudotime, reference = gene1, target = gene2)
} 