#' Plot heatmap of marker genes
#'
#' @description Plot heatmap of gene expression for single cell data
#'
#' @param object The \code{SingleCellExperiment} S4 object
#' @param markers the data.frame, list or character vector of genes/rows to be plotted
#' @param coldata.to.include the colData columns to be included in the col_anno of the heatmap
#' @param cluster_rows cluster_rows
#' @param cluster_cols cluster_cols
#' @param show_rownames show_rownames
#'
#' @importFrom pheatmap pheatmap
#' @import SingleCellExperiment
#'
#' @export
PlotExprHeatmap <- function(object, markers = NULL, coldata.to.include = NULL, cluster_rows = F, cluster_cols = F, show_rownames = T, use.gap = T, ...){

  if(mode(markers) %in% c("list", "data.frame") & any(colnames(markers) == "gene")){ # input marker genes is a data frame
    markers <- markers[markers$gene %in% rownames(object),]
    plot.data <- logcounts(object)[match(markers$gene, rownames(object)),]
  } else if (mode(markers) == "character") { # input marker genes is a vector
    markers <- markers[markers %in% rownames(object)]
    plot.data <- logcounts(object)[match( markers, rownames(object)),]
  }

  if(!is.null(coldata.to.include)){
    coldata.to.include <- as.character(coldata.to.include[coldata.to.include %in% colnames(colData(object))])
    coldata.to.include <- c("final.clusters",coldata.to.include)
  } else {
    coldata.to.include <- "final.clusters"
  }

  if(length(coldata.to.include) > 1){
    colanno <- data.frame (colData(object)[,match(coldata.to.include,colnames(colData(object)) )])
    colnames(colanno)[1] <- "clusters"
    colanno$clusters <- factor(colanno$clusters)
  } else {
    colanno <- data.frame (clusters = colData(object)$final.clusters)
  }
  rownames(colanno) <- colnames(object)

  # colanno <- colanno[order(colanno$order),, drop = F]
  # colanno$clusters <- as.numeric(as.character(colanno$clusters))
  colanno <- colanno[order(colanno$clusters, decreasing = F),]
  # colanno <- colanno[,1:2]
  colanno$clusters <- factor(colanno$clusters, levels = unique(colanno$clusters))

  plot.data <- plot.data[,match(rownames(colanno), colnames(plot.data))]

  if(use.gap == TRUE){
    gap <- which(c(colanno$clusters, 99)!= c(99, colanno$clusters))
    gap <- gap[gap <= nrow(colanno) & gap > 1]
    gap <- gap - 1
  } else {
    gap <- 0
  }

  rowanno  <- NULL
  if(mode(markers) %in% c("list", "data.frame") & any(colnames(markers) == "cluster")){ # input is a data frame
    rowanno <- data.frame (clusters = markers$cluster)
    rownames(rowanno) <-  markers$gene
    plot.data <- plot.data[match(rownames(rowanno), rownames(plot.data)),]
    rowanno$clusters <- as.character(rowanno$clusters)
    # rowanno <- rowanno[,1, drop = F]
  }

  plot.data_scale <- plot.data

  # plot.data_scale <- scale(plot.data_scale, center = T, scale = T)
  r_mean <- apply(plot.data_scale, 1, mean)

  for(i in 1:nrow(plot.data_scale)){ # scale data
    plot.data_scale[i,] <- plot.data[i,] - r_mean[i]
  }
  r_max <- apply(abs(plot.data_scale), 1, max)
  for(i in 1:nrow(plot.data_scale)){ # scale data
    plot.data_scale[i,] <- plot.data_scale[i,]/r_max[i] * 5
  }

  print(dim(plot.data_scale))
  print(colanno[1:5,])
  my_colours <- colorRampPalette(c("steelblue3", "white", "firebrick3"))(200)
  if(is.null(rowanno)) {
    pheatmap(plot.data_scale, color = my_colours,
             annotation_col = colanno, scale = "none", gaps_col = gap,
             cluster_rows = cluster_rows, cluster_cols = cluster_cols,
             show_colnames = F, show_rownames = show_rownames, ...)
  } else {
    pheatmap(plot.data_scale, color = my_colours,
             annotation_row = rowanno, annotation_col = colanno, scale = "none", gaps_col = gap,
             cluster_rows = cluster_rows, cluster_cols = cluster_cols,
             show_colnames = F, show_rownames = show_rownames, ...)
  }


}
