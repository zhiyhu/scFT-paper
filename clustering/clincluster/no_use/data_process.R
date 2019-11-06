#' Prepare data before the clustering
#'
#' @description Prepare data before the clustering, including normalisation, identification of high variance genes and regressing unwanted variables
#'
#' @param object an object of 'SingleCellExperiment' class
#' @param col.for.cluster the colnames that contains the patient/batch information for initial clustering
#### @param do.regress whether to regress out the effect of sequencing depth on data
#' @param log.scale.factor log.scale.factor
#' @param verbose verbose
#' @param do.scale do.scale
#' @param do.centre do.centre
#'
#' @return object an object of 'SingleCellExperiment' class
#'
#### @importFrom Seurat LogNormalize
#' @import SingleCellExperiment
#' @importFrom scater calculateCPM
#'
#' @export
PrepareData <- function(object, col.for.cluster,
                        # do.regress = F,
                        log.scale.factor = 100000,
                        verbose = T,
                        do.scale = F,
                        do.centre = T)
{
  cell.info <- object@colData
  expr <- counts(object)
  if(all(rownames(cell.info) == colnames(expr))){

    if(col.for.cluster %in% colnames(object@colData)){
      tmp <- cell.info[,which(colnames(object@colData) == col.for.cluster)]
      object@colData$group.for.cluster <- as.character(tmp) ## extract the information of batches
    } else {
      print("col.for.cluster does not exist. Please check")
    }

    # lognormalisation
    if(verbose){
      print("Normalising the data...")
    }
    # norm.expr <- LogNormalize(data = expr, scale.factor = log.scale.factor, verbose = verbose)
    norm.expr <- log1p(calculateCPM(expr)*log.scale.factor/1e6) # log(CPM+1)
    logcounts(object) <- as.matrix(norm.expr)

    logcounts(object) <- as.matrix(norm.expr)
    gc()

    # if(do.regress){ logcounts(object) <- RegressOutUnwanted.seurat(object, verbose = verbose) }
    if(verbose){
      print("Centre the data")
    }
    scalecounts <- scale(logcounts(object), scale = do.scale, center = do.centre)
    normcounts(object) <- scalecounts
    return(object)

  }else{
    print("rownames(cell) do not match colnames(expr). Please check.")
  }
}

#' Identify high variance genes (adopted from Seurat)
#'
#' @description Identify high variance genes from normalised log count matrix
#'
#' @param object an object of 'SingleCellExperiment' class
#' @param verbose show the progress bar or not (default: \code{TRUE})
#' @param genes.use genes.use
#' @param num.bin num.bin
#' @param mean.low.cutoff mean.low.cutoff
#' @param mean.high.cutoff mean.high.cutoff
#' @param dispersion.low.cutoff dispersion.low.cutoff
#' @param dispersion.high.cutoff dispersion.high.cutoff
#'
#' @import ggplot2
#'
#' @return A charactor vector that contains the high variance genes
#' @references Seurat
#'
#' @export
HighVarGenes <- function(object,
                         verbose = F,
                         genes.use = NULL,
                         num.bin = 20,
                         mean.low.cutoff = 0.1, mean.high.cutoff = 8,
                         dispersion.low.cutoff = 1, dispersion.high.cutoff = Inf) {

  data <- logcounts(object = object) # log expression
  # identify high variance genes
  if(is.null(genes.use)) {
    genes.use <- rownames(data) # the genes of rownames of data
  } else {
    genes.use <- genes.use[genes.use %in% rownames(data)]
    data <- data[match(genes.use, rownames(data)), ]
  }
  # set the mean and dispersion of genes to zero.
  gn.mean <- rep(0, length(genes.use))
  gn.dispr <- rep(0, length(genes.use))
  gn.dispr.scaled <- rep(0, length(genes.use))
  names(gn.mean) <- names(gn.dispr) <- names(gn.dispr.scaled) <- genes.use

  if(verbose){
    print("Finding high variance genes...")
  }

  gn.mean <- apply(data, MARGIN = 1, FUN = ExprMean) # mean of all genes across cells
  gn.dispr <- apply(data, MARGIN = 1, FUN = ExprVar) # dispersion of all genes

  gn.dispr[is.na(gn.dispr)] <- 0 # set NA to 0
  gn.mean[is.na(gn.mean)] <- 0 # set NA to 0

  index_bin <- cut(gn.mean, breaks = num.bin)
  names(index_bin) <- names(gn.mean)

  mean_bin <- tapply(gn.dispr, INDEX = index_bin, FUN = mean) # calculate mean for dispersion of each bin
  sd_bin <- tapply(gn.dispr, INDEX = index_bin, FUN = sd) # calculate SD for dispersion of each bin

  # scale dispersion of each bin by the mean and dispersion
  gn.dispr.scaled <- (gn.dispr - mean_bin[as.numeric(index_bin)]) / sd_bin[as.numeric(index_bin)]
  gn.dispr.scaled[is.na(gn.dispr.scaled)] <- 0 # set NA to 0
  names(gn.mean) <- names(gn.dispr) <- names(gn.dispr.scaled) <- rownames(data)

  var.df <- data.frame(gn.mean = gn.mean,
                       gn.dispr = gn.dispr,
                       gn.dispr.scaled = gn.dispr.scaled)
  rownames(var.df) <- rownames(data)

  var.genes <- names(gn.mean)[which((gn.mean > mean.low.cutoff) &
                                      (gn.mean < mean.high.cutoff) &
                                      (gn.dispr.scaled > dispersion.low.cutoff) &
                                      (gn.dispr.scaled < dispersion.high.cutoff))]
  if(verbose){
    print(paste(length(var.genes),"high variance genes are found.", sep = " "))
  }
  rowData(object)$high.var <- F
  rowData(object)$high.var[rownames(object) %in% var.genes] <- T
  rowData(object)$gene.mean <- gn.mean
  rowData(object)$gene.dispersion <- gn.dispr
  rowData(object)$gene.dispersion.scaled <- gn.dispr.scaled
  return(object)
}

# Calculate the mean of log-expression
ExprMean <- function(x) {
    return(log1p(mean(expm1(x))))
}

# Calculate the ratio of variance/mean of log-expression
ExprVar <- function(x) {
    return(log(var(expm1(x)) / mean(expm1(x))))
}


#' Prepare data before the clustering
#'
#' @description Prepare data before the clustering, including normalisation, identification of high variance genes and regressing unwanted variables
#'
#' @param object an object of 'SingleCellExperiment' class
#' @param pca.max.component pca.max.component
#' @param use.hvg use high variance genes
#' @return object an object of 'SingleCellExperiment' class
#'
#' @importFrom stats prcomp
#' @import SingleCellExperiment
#' @export
calculatePCA <- function(object, use.hvg = T,
                         pca.max.component = 50){

  scalecounts <- normcounts(object)
  # calculate PCA from high variance genes
  if(use.hvg == T) {
    pca <- prcomp(t(scalecounts[rowData(object)$high.var,]), rank. = pca.max.component)
  } else {
    pca <- prcomp(t(scalecounts), rank. = pca.max.component)
  }

  ttl.var <- sum(pca$sdev ^ 2) # total variance
  pct.Var <- pca$sdev ^ 2 / ttl.var # percentage variance
  pcs <- pca$x
  attr(pcs, "percentVar") <- pct.Var

  reducedDim(object, "PCA") <- pcs

  return(object)
}

