#' Plot heatmap of the gene expression
#'
#' XXXXXXX
#'
#' @name PlotHeatmap
#'
#' @param object an object of 'SCESet' class
#' @param show_pdata a vector of colnames of the pData(object) table. Default is NULL.
#' If not NULL will add pData annotations to the columns of the output matrix
#'
#' @importFrom pheatmap pheatmap
#'
#' @export
PlotHeatmap <- function(dataset, show_pdata = NULL) {
  if (is.null(dataset$cell.info$final.clusters)) {
    warning(paste0("Please run clustering first!"))
  }
  else{
    add_ann_col <- FALSE
    ann <- NULL
    if (!is.null(show_pdata)) {
      ann <- make_col_ann_for_heatmaps(object, show_pdata)
      if (!is.null(ann)) {
        add_ann_col <- TRUE
        # make same names for the annotation table
        rownames(ann) <- colnames(consensus)
      }
    }
    do.call(pheatmap::pheatmap, c(list(consensus, cluster_rows = hc, cluster_cols = hc, cutree_rows = k,
                                       cutree_cols = k, show_rownames = FALSE, show_colnames = FALSE), list(annotation_col = ann)[add_ann_col]))

  }
}


#' Plot hierarchical tree of the final clustering
#'
#' XXXXXXX
#'
#'
#' @param dataset dataset
#' @param col.border col.border
#'
#' @export
PlotHierarchicalTree <- function(dataset, col.border = "red")
{
  plot(dataset[[6]])
  rect.hclust(dataset[[6]], k = dataset[[7]], border = col.border)
}
