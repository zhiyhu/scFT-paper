#' Initially cluster within groups
#'
#' @description
#'
#' @param object the object from \code{PrepareData}
#' @param method the method used for the initial clustering. The available parameters are "spectral" (default and recommended), "ap" (affinity propagation) and "kmeans" .
#' @param k number of clusters if \code{method = "kmeans"} or \code{method = "spectral"}.
#' @param spec.method method for spectral clustering
#' @param n.neighbor number of knn neighbours to be used in the kknn method
#' @param display.progress display the progress bar (default = T)
#' @param kernel the kernel function used in computing the affinity matrix. The default is the Gaussian kernal \code{"rbfdot"}.
#' @param kpar a character string or the list of hyper-parameters. The default is  \code{"automatic"}.
#' @param nystrom.red use nystrom method to calculate eigenvectors. Recommanded for large number of cells. (default: \code{FALSE})
#' @param nystrom.sample number of data points to use for estimating the eigenvalues when using the nystrom method. (default : \code{dim(x)[1]/6})
#' @param mod.sample proportion of data to use when estimating sigma (default: 0.75).
#' @param iterations the maximum number of iterations allowed.
#' @param use.pca use PCA for clustering or not
#' @param ncomponents n PCA
#'
#' @details
#'
#' @import SingleCellExperiment
#' @importFrom scater runPCA
#'
#' @return A \code{SingleCellExperiment} S4 object
#'
#' @export


# InitialCluster <- function(object, method = "spectral", display.progress = T,
#                            ncomponents = 10,
#                            k = 3, regressed = F,
#                            use.var.genes = T,
#                            use.pca = T,
#                            n.pca = 7,
#                            spec.method = "kknn",
#                            n.neighbor = 5,
#                            # parameters for spectral clustering
#                            kernel = "rbfdot", kpar = "automatic",
#                            nystrom.red = FALSE, nystrom.sample = NULL,
#                            mod.sample = 0.75, iterations = 400,...){
#     
#     # if(use.pca == F){ # not using pca for clustering
#     #     scalecounts <-  as.matrix(logcounts(object))
#     #     scalecounts <- scale(logcounts, scale = do.scale, center = do.centre)# scale on the genes
#     #     if(use.var.genes){ # use high varia
#     #         scalecounts <- scalecounts[rowData(object)$high.var == T,] # matrix for initial clustering
#     #     }
#     # } else if(use.pca == T){ # using pca for clustering
#     #     scalecounts <- t(object@reducedDims$PCA[,1:n.pca])
#     # }
#     groups <- names(table(object@colData$group.for.cluster))
#     n.group <- length(groups) # number of patients
#     
#     for(i in 1:n.group)
#     {
#         #extract the parts for clustering
#         subset <- object@colData$group.for.cluster %in% groups[i]
#         subobj <- object[,subset]
#         
#         cell <- subobj@colData # need to transfer it into a data frame
#         cell$Sample <- rownames(cell)
#         
#         # Run PCA reduction
#         sceset <- runPCA(object = subobj, ncomponents = 20, exprs_values = "logcounts", rand_seed = 12345,
#                          feature_set = rownames(subobj)[rowData(subobj)$high.var == TRUE])
#         
#         plot(1:50, (attr(subobj@reducedDims$PCA, "percentVar")[1:50])*100, pch = 20, xlab = "PC", ylab = "Standard Deviation of PC", main = groups[i])
#         
#         scalecounts <- t(subobj@reducedDims$PCA[,1:n.pca])
#         
#         cell <- InitialCluster_spectral(scalecounts = scalecounts, cell = cell, k = k[i],
#                                         do.scale = do.scale, do.centre = do.centre,
#                                         method = spec.method, nn = n.neighbor,
#                                         kernel = kernel, kpar = kpar, nystrom.red = nystrom.red, nystrom.sample = nystrom.sample,
#                                         mod.sample = mod.sample, iterations = iterations,
#                                         display.progress = display.progress)
#         
#         colData(object)$initial.cluster[match(rownames(cell), rownames(object@colData))] <- cell$initial.cluster
#         
#     }
#     
#     
#     return(object)
# }
# 


InitialCluster <- function(object, method = "spectral", display.progress = T,
                           ncomponents = 10,
                           k = 3, regressed = F,
                           use.var.genes = T,
                           use.pca = T,
                           n.pca = 7,
                           spec.method = "kknn",
                           n.neighbor = 5,
                           # parameters for spectral clustering
                           kernel = "rbfdot", kpar = "automatic",
                           nystrom.red = FALSE, nystrom.sample = NULL,
                           mod.sample = 0.75, iterations = 400,...){

  cell <- object@colData # need to transfer it into a data frame
  cell$Sample <- rownames(cell)

  if(use.pca == F){ # not using pca for clustering
    scalecounts <-  as.matrix(logcounts(object))
    scalecounts <- scale(logcounts, scale = do.scale, center = do.centre)# scale on the genes
    if(use.var.genes){ # use high varia
      scalecounts <- scalecounts[rowData(object)$high.var == T,] # matrix for initial clustering
    }
  } else if(use.pca == T){ # using pca for clustering
     scalecounts <- t(object@reducedDims$PCA[,1:n.pca])
  }

  print("Got data ready")
  if(method == "kmeans"){
    cell <- InitialCluster_kmeans(scalecounts, cell, k)
  } else if(method == "spectral") {

    cell <- InitialCluster_spectral(scalecounts = scalecounts, cell = cell, k,do.scale = do.scale, do.centre = do.centre,
                                    method = spec.method, nn = n.neighbor,
                                    kernel = kernel, kpar = kpar, nystrom.red = nystrom.red, nystrom.sample = nystrom.sample,
                                    mod.sample = mod.sample, iterations = iterations,
                                    display.progress = display.progress)
  } else if(method == "ap") {
    cell <- InitialCluster_AP(scalecounts, cell)
  } else {
    print("Unrecognised method. Default method (spectral clustering) is used.")
    cell <- InitialCluster_spectral(scalecounts, cell, k)
  }
  colData(object)$initial.cluster <- cell$initial.cluster

  return(object)
}



#' Calculate distance matrix by differential expression analysis
#'
#' @description Calculate distance matrix by differential expression (DE) analysis
#'
#' @param object The \code{SingleCellExperiment} S4 object returned from \code{InitialCluster()}
#' @param use.var.genes use high variance genes only in DE analysis (default: \code{TRUE}). Speed up the calculation and more stable.
#' @param p.value The p.value to screen the differentially expressed genes and to be passed to the \code{limma::topTable} （default: 0.01)
#' @param n.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores. The default is generated by \code{detectCores()}.
#' @param low.pct cutoff needed to identify DE genes. The lower cutoff of percentage of cells expressing the genes.
#' @param log.FC the least log-fold change needed for DE genes.
#' @param display.progress show the progress bar or not (default: \code{TRUE})
#' @param sampling sampling or not
#'
#' @details Calculate distance matrix by differential expression analysis
#'
#' @return return the distance matrix for final clustering. Use it as the input of \code{FinalCluster()}
#'
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom parallel mclapply
#' @import SingleCellExperiment
#'
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), pp. e47.
#'
#' @seealso limma voom
#'
#' @export
CalculateDEMatrix <- function(object,
                              use.var.genes = T,
                              p.value = 0.05,
                              low.pct = 0.3, # expr.quantile = 0.25,
                              log.FC = 0.6,
                              sampling = F,
                              n.sampling = 1,
                              display.progress = T) {
  cell <- colData(object)
  initial.cluster <- cell$initial.cluster

  initial.cluster.names <- names(table(initial.cluster)) # the names of initial clusters
  # calculate the DE distance matrix
  m_de <- matrix(0, nrow = length(initial.cluster.names), # DE matrix
                 ncol = length(initial.cluster.names))
  colnames(m_de) <- rownames(m_de) <- initial.cluster.names
  nc <- ncol(m_de)

  m_de_df <- data.frame(i = rep(1:nc, nc), j = rep(1:nc,each = nc)) # dataframe for cycling
  m_de_df <- m_de_df[m_de_df$j > m_de_df$i,]
  rownames(m_de_df) <- 1:nrow(m_de_df)

  m_de_list <- as.list(as.data.frame(t(as.matrix(m_de_df))))
  res_all <- list()

  # if(is.null(n.cores)){ n.cores <- detectCores()}
  if(display.progress ==T & sampling == T){
    print("Calculate distance matrix")
    pb <- txtProgressBar(min = 0, max = length(m_de_list)*n.sampling, style = 3)
  } else if((display.progress ==T & sampling == F)){
    print("Calculate distance matrix")
    pb <- txtProgressBar(min = 0, max = length(m_de_list), style = 3)
  }

  for (k in 1:length(m_de_list)) {
    i <- initial.cluster.names[m_de_list[[k]][1]] # cluster i
    j <- initial.cluster.names[m_de_list[[k]][2]] # cluster j

    if (sampling == T) { # sampling 80% of the samples

      idx_i <- which(colData(object)$initial.cluster %in% c(i))
      idx_j <- which(colData(object)$initial.cluster %in% c(j))
      lgt <- min(floor(length(idx_i) * 0.8), floor(length(idx_j) * 0.8)) #size of each group
      res.n <- c()

      for (ns in 1:n.sampling) { # do sampling DE analysis for n.sampling times
        idx1 <- sample(x = idx_i, size = lgt )
        idx2 <- sample(x = idx_j, size = lgt )
        index <- c(idx1, idx2)

        if(use.var.genes){ # matrix is the expr logcounts
          matrix <- expm1(logcounts(object))[rowData(object)$high.var == TRUE, index] # expr matrix
        } else {
          matrix <- expm1(logcounts(object))[ ,index]
        }
        sample.info <- colData(object)$initial.cluster[index]
        dim(matrix)

        dge <- edgeR::DGEList(counts = matrix) # make a edgeR object
        design <- model.matrix(~sample.info)
        colnames(design) <- names(table(sample.info))

        v <- limma::voom(dge, design, plot=F)
        fit <- limma::lmFit(v, design) # Fit linear model for each gene given a series of arrays.
        fit <- limma::eBayes(fit)
        index2.1 <- rowSums(matrix[,sample.info %in% unique(sample.info)[1]] > 0)/ncol(matrix[,sample.info %in% unique(sample.info)[1]])
        index2.2 <- rowSums(matrix[,sample.info %in% unique(sample.info)[2]] > 0)/ncol(matrix[,sample.info %in% unique(sample.info)[2]])
        index2 <- rowMax(cbind(index2.1, index2.2))
        names(index2) <- rownames(matrix)

        # res.all <-  limma::topTable(fit, coef = 2, number = Inf)
        # res.all$ratio <- index2[match(rownames(res.all), names(index2))]
        res <- limma::topTable(fit, coef = 2, number = Inf, p.value = p.value)
        res$ratio <- index2[match(rownames(res), names(index2))]
        res$ratio1 <- index2.1[match(rownames(res), names(index2.1))]
        res$ratio2 <- index2.2[match(rownames(res), names(index2.2))]

        res <- res[abs(as.numeric(res$logFC)) > log.FC,]
        res <- res[res$ratio > low.pct,]
        res.n <- c(res.n, nrow(res))

        if(display.progress){
          setTxtProgressBar(pb, n.sampling*k+ns)
        }
      }

      m_de_df$res[k] <- mean(res.n)

    } else if (sampling == F) { # do not use sampling to choose the DE groups

      index <- which(colData(object)$initial.cluster %in% c(i,j))
      if(use.var.genes){ # matrix is the expr logcounts
        matrix <- expm1(logcounts(object))[rowData(object)$high.var == TRUE, index] # expr matrix
      } else {
        matrix <- expm1(logcounts(object))[ ,index]
      }
      sample.info <- colData(object)$initial.cluster[index]

      dge <- edgeR::DGEList(counts = matrix) # make a edgeR object
      design <- model.matrix(~sample.info)
      colnames(design) <- names(table(sample.info))

      v <- limma::voom(dge, design, plot=F)
      fit <- limma::lmFit(v, design) # Fit linear model for each gene given a series of arrays.
      fit <- limma::eBayes(fit)
      index2.1 <- rowSums(matrix[,sample.info %in% unique(sample.info)[1]] > 0)/ncol(matrix[,sample.info %in% unique(sample.info)[1]])
      index2.2 <- rowSums(matrix[,sample.info %in% unique(sample.info)[2]] > 0)/ncol(matrix[,sample.info %in% unique(sample.info)[2]])
      index2 <- rowMax(cbind(index2.1, index2.2))
      names(index2) <- rownames(matrix)

      res.all <-  limma::topTable(fit, coef = 2, number = Inf)
      res.all$ratio <- index2[match(rownames(res.all), names(index2))]

      res <- limma::topTable(fit, coef = 2, number = Inf, p.value = p.value)
      res$ratio <- index2[match(rownames(res), names(index2))]
      res$ratio1 <- index2.1[match(rownames(res), names(index2.1))]
      res$ratio2 <- index2.2[match(rownames(res), names(index2.2))]

      res_all[[i]] <- res

      res <- res[abs(as.numeric(res$logFC)) > log.FC,]
      res <- res[res$ratio > low.pct,]

      if(length(nrow(res))==0)
      {
        m_de_df$res[k] <- 0
      }else{
        m_de_df$res[k] <- nrow(res)
      }

      if(display.progress){
        setTxtProgressBar(pb, k)
      }

     }
  }
  # print("[2/3] Generating distance matrix...")
  # res <- parallel::mclapply(m_de_list, RunLimmaVoom, matrix = matrix, p.value = p.value,
  # length.contrast = length(initial.cluster.names), mc.cores = n.cores)
  # res <- unlist(res); m_de_df$res <- as.numeric(res)

  # print("[3/3] Finishing...")
  for(k in 1:nrow(m_de_df)){ # transfer DE matrix from data.frame to a matrix
    m_de[m_de_df$i[k], m_de_df$j[k]] <- as.numeric(m_de_df$res[k])
  }
  DE.matrix <- as.dist(as.matrix(m_de)+t(as.matrix(m_de))) # transfer to dist object

  if(display.progress){
    close(pb)
  }

  return(list(DE.matrix = DE.matrix,
              # DE.res = res_all,
              comparison = m_de_list))
}


#' Final Clustering of the initial clusters
#'
#' @description Run the clustering of the initial clusters based on the pre-calculated distance matrix
#'
#' @param object The \code{SingleCellExperiment} S4 object returned from \code{InitialCluster()}
#' @param DE.matrix The returned value from \code{CalculateDEMatrix()}
#' @param k.cluster The number of final clusters
#' @param plot Whether to plot the hierarchical tree of the final clustering (default: TRUE)
#'
#' @details Run the clustering of the initial clusters based on the pre-calculated distance matrix
#'
#' @return A \code{SingleCellExperiment} S4 object that contains that final clustering results
#'
#' @examples
#'
#' @export
FinalCluster <- function(object, DE.matrix, k.cluster, plot = TRUE){

  hc <- hclust(as.dist(DE.matrix))
  hc.cluster <- cutree(hc, k = k.cluster)

  colData(object)$final.clusters <- hc.cluster[match(colData(object)$initial.cluster, names(hc.cluster))]
  colData(object)$final.clusters <- as.factor(colData(object)$final.clusters)

  if(plot == T){
    plot(hc)
    rect.hclust(hc, k = k.cluster, border = "red")
  }
  return(object)
}

#' Identify marker genes
#'
#' @param object The \code{SingleCellExperiment} S4 object
#' @export
#'
GetMarkers <- function(object)
{
  cell <- dataset[[3]]
  dge <- dataset[[4]]

  final.cluster <- cell$final.clusters

  design <- model.matrix(~final.cluster)
  colnames(design) <- names(table(final.cluster))
  print("[1/3] Fitting linear model...")

  v <- limma::voom(dge, design, plot=F)
  fit <- limma::lmFit(v, design) # Fit linear model for each gene given a series of arrays.
  # return MArrayLM fitted model object
  initial.cluster.names <- names(table(initial.cluster)) # the names of initial clusters

  return(object)

}

#' Calculate the DE matrix between patients
#'
#' @description Calculate distance matrix by differential expression (DE) analysis
#'
#' @param object The \code{SingleCellExperiment} S4 object returned from \code{InitialCluster()}
#' @param use.var.genes use high variance genes only in DE analysis (default: \code{TRUE}). Speed up the calculation and more stable.
#' @param p.value The p.value to screen the differentially expressed genes and to be passed to the \code{limma::topTable} （default: 0.01)
#' @param n.cores the number of cores to use, i.e. at most how many child processes will be run simultaneously. Must be at least one, and parallelization requires at least two cores. The default is generated by \code{detectCores()}.
#' @param low.pct cutoff needed to identify DE genes. The lower cutoff of percentage of cells expressing the genes.
#' @param log.FC the least log-fold change needed for DE genes.
#' @param display.progress show the progress bar or not (default: \code{TRUE})
#'
#' @details Calculate distance matrix by differential expression analysis
#'
#' @return return the distance matrix for final clustering. Use it as the input of \code{FinalCluster()}
#'
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#' @importFrom parallel mclapply
#'
#' @references Ritchie ME, Phipson B, Wu D, Hu Y, Law CW, Shi W and Smyth GK (2015). “limma powers differential expression analyses for RNA-sequencing and microarray studies.” Nucleic Acids Research, 43(7), pp. e47.
#'
#' @seealso limma voom
#' @export
CalculateBgDEMatrix <- function(object,
                                use.var.genes = T,
                                p.value = 0.05,
                                low.pct = 0.3, # expr.quantile = 0.25,
                                log.FC = 0.8,
                                display.progress = T) {
  cell <- colData(object)
  initial.cluster <- cell$group.for.cluster

  initial.cluster.names <- names(table(initial.cluster)) # the names of initial clusters
  # calculate the DE distance matrix
  m_de <- matrix(0, nrow = length(initial.cluster.names), # DE matrix
                 ncol = length(initial.cluster.names))

  colnames(m_de) <- rownames(m_de) <- initial.cluster.names
  nc <- ncol(m_de)

  m_de_df <- data.frame(i = rep(1:nc, nc), j = rep(1:nc,each = nc)) # dataframe for cycling
  m_de_df <- m_de_df[m_de_df$j > m_de_df$i,]
  rownames(m_de_df) <- 1:nrow(m_de_df)

  m_de_list <- as.list(as.data.frame(t(as.matrix(m_de_df))))
  res_all <- list()

  # if(is.null(n.cores)){ n.cores <- detectCores()}
  if(display.progress){
    print("Calculate distance matrix")
    pb <- txtProgressBar(min = 0, max = length(m_de_list), style = 3)
  }

  for (k in 1:length(m_de_list)) {
    i <- initial.cluster.names[m_de_list[[k]][1]] # cluster i
    j <- initial.cluster.names[m_de_list[[k]][2]] # cluster j

    index <- which(colData(object)$group.for.cluster %in% c(i,j))

    if(use.var.genes){ # matrix is the expr logcounts
      matrix <- expm1(logcounts(object))[rowData(object)$high.var == TRUE, index] # expr matrix
    } else {
      matrix <- expm1(logcounts(object))[ ,index]
    }

    sample.info <- colData(object)$group.for.cluster[colData(object)$group.for.cluster %in% c(i,j)]

    dge <- edgeR::DGEList(counts = matrix) # make a edgeR object

    design <- model.matrix(~sample.info)
    colnames(design) <- names(table(sample.info))

    v <- limma::voom(dge, design, plot=F)
    fit <- limma::lmFit(v, design) # Fit linear model for each gene given a series of arrays.
    fit <- limma::eBayes(fit)
    index2.1 <- rowSums(matrix[,sample.info %in% unique(sample.info)[1]] > 0)/ncol(matrix[,sample.info %in% unique(sample.info)[1]])
    index2.2 <- rowSums(matrix[,sample.info %in% unique(sample.info)[2]] > 0)/ncol(matrix[,sample.info %in% unique(sample.info)[2]])
    index2 <- rowMax(cbind(index2.1, index2.2))
    names(index2) <- rownames(matrix)

    res.all <-  limma::topTable(fit, coef = 2, number = Inf)
    res.all$ratio <- index2[match(rownames(res.all), names(index2))]

    res <- limma::topTable(fit, coef = 2, number = Inf, p.value = p.value)
    res$ratio <- index2[match(rownames(res), names(index2))]
    res$ratio1 <- index2.1[match(rownames(res), names(index2.1))]
    res$ratio2 <- index2.2[match(rownames(res), names(index2.2))]

    res_all[[i]] <- res

    res <- res[abs(as.numeric(res$logFC)) > log.FC,]
    res <- res[res$ratio > low.pct,]

    if(length(nrow(res))==0)
    {
      m_de_df$res[k] <- 0
    }else{
      m_de_df$res[k] <- nrow(res)
    }
    if(display.progress){
      setTxtProgressBar(pb, k)
    }
  }

  for(k in 1:nrow(m_de_df)){
    m_de[m_de_df$i[k], m_de_df$j[k]] <- as.numeric(m_de_df$res[k])
  }
  DE.matrix <- as.dist(as.matrix(m_de)+t(as.matrix(m_de)))

  if(display.progress){
    close(pb)
  }

  return(list(DE.matrix = DE.matrix, DE.res = res_all, comparison = m_de_list))
}

