#' Identify the differentially expressed genes between two populations
#'
#' @description Identify the differentially expressed genes between two populations
#'
#' @param object The \code{SingleCellExperiment} S4 object
#' @param log.FC the cutoff of log-fold change needed for DE genes.
#' @param display.process display the processing bar
#' @param strict whether the marker should be the marker of the only cluster against all the other clusters
#'
#' @details Identify the marker genes for each final cluster by limma voom
#'
#' @return data frame of marker genes and results of differential expression analysis
#'
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @import SingleCellExperiment
#'
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), pp. e47.
#'
#' @seealso IdentifyDEGenes
#'
#' @export
IdentifyMarkers <- function(object, log.FC = 0.4, pos.only = T, display.progress = T, strict = T)
{
  if(is.null(object$final.clusters)){
    print("Please run the clustering steps first")
    return(object)
  }

  names.clusters <- unique(object$final.clusters)
  names.clusters <- names.clusters[order(names.clusters, decreasing = F)]
  marker_df <- NULL

  if(strict == T){ # strict is T
    for(i in 1:length(names.clusters))
    {
      marker_tmp <- c()
      cluster <- names.clusters[i]
      other.clusters <- names.clusters[-i]
      if(display.progress){
        print(paste("Cluster",cluster,"Makers"))
        pb <- txtProgressBar(min = 0, max =  (length(names.clusters)), style = 3)
      }

      for(j in 1:length(other.clusters)) # do DE analyss between one cluster and every other cluster
      {
        de_tmp <- IdentifyDEGenes(object, print.process = F, # DE gene list
                                  idx.2 = which(object$final.clusters == cluster),
                                  idx.1 = which(object$final.clusters == other.clusters[j]),
                                  log.FC = log.FC)

        de_tmp$gene <- rownames(de_tmp)
        de_tmp <- de_tmp[de_tmp$logFC > 0,] # subset by upregulated genes in designated cluster

        # merge the data with previous data
        if(length(marker_tmp) == 0){
          marker_tmp <- de_tmp$gene
        } else {
          marker_tmp <- marker_tmp[marker_tmp %in% de_tmp$gene]
        }
        if(display.progress){
          setTxtProgressBar(pb, j)
        }
      }

      # identify the DE genes between the backgroup clusters and the designated cluster
      de <- IdentifyDEGenes(object, print.process = F,
                            idx.1 = which(object$final.clusters !=cluster),
                            idx.2 = which(object$final.clusters == cluster),
                            log.FC = log.FC)
      if(display.progress){
        setTxtProgressBar(pb, j+1)
      }
      de <- de[order(de$logFC, decreasing = T),] # order the table by logFC
      de$gene <- rownames(de)
      de$cluster <- cluster
      de <- de[de$gene %in% marker_tmp,] # choose the markers

      # merge the data with previous data
      if(is.null(marker_df)){
        marker_df <- de
      } else {
        marker_df <- rbind(marker_df, de)
      }
      if(display.progress){
        close(pb)
      }
    }
  } else { # if strict is F
    for(i in 1:length(names.clusters)) # for every cluster
    {
      marker_tmp <- c()
      cluster <- names.clusters[i]

      if(display.progress){
        print(paste("Cluster",cluster,"Makers"))
      }
      # identify the DE genes between the backgroup clusters and the designated cluster
      de <- IdentifyDEGenes(object, print.process = display.progress,
                            idx.1 = which(object$final.clusters !=cluster),
                            idx.2 = which(object$final.clusters == cluster),
                            log.FC = log.FC)

      de <- de[order(de$logFC, decreasing = T),] # order the table by logFC
      de$gene <- rownames(de)
      de$cluster <- cluster

      # merge the data with previous data
      if(is.null(marker_df)){
        marker_df <- de
      } else {
        marker_df <- rbind(marker_df, de)
      }
    }
  }


  return(marker_df)
}
