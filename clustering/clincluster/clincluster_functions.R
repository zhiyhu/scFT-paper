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
                        do.centre = T) {
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

#' Initially cluster within groups by spectral clustering
#'
#' @param scalecounts the count matrix after centering and scaling by \code{scale()}
#' @param cell a data.frame that contains cell information; its rownames must be same with the colnames of scalecounts
#' @param k number of clusters
#' @param display.progress display the progress bar (default = T)
#' @param kernel the kernel function used in computing the affinity matrix. The default is the Gaussian kernal \code{"rbfdot"}.
#' @param kpar a character string or the list of hyper-parameters. The default is  \code{"automatic"}.
#' @param nystrom.red use nystrom method to calculate eigenvectors. Recommanded for large number of cells. (default: \code{FALSE})
#' @param nystrom.sample number of data points to use for estimating the eigenvalues when using the nystrom method. (default : \code{dim(x)[1]/6})
#' @param mod.sample proportion of data to use when estimating sigma (default: 0.75).
#' @param iterations the maximum number of iterations allowed.
#' @param method kernlab or kknn
#'
#' @importFrom kernlab specc
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom kknn specClust
#'
#' @references 12345
#' @export
#' @seealso kernlab::specc
#' 
InitialCluster_spectral <- function(scalecounts, cell, k, display.progress = T,
                                    kernel = "rbfdot", kpar = "automatic",
                                    nystrom.red = FALSE, nystrom.sample = NULL,
                                    mod.sample = 0.85, iterations = 200,
                                    method = "kernlab",
                                    ...){
    groups <- names(table(cell$group.for.cluster))
    n.group <- length(groups) # number of patients
    cell$initial.cluster <- NA
    
    if(display.progress){
        print("Initial clustering")
        pb <- txtProgressBar(min = 0, max = n.group, style = 3)# create progress bar
    }
    
    if (length(k) != 1 & length(k) != n.group & (!is.null(k))) { # wrong length of k
        
        print("Incorrect length of k. Please check...")
        return(cell)
        
    } else if (length(k) == 1 | is.null(k)) { # only a uniform value for number of clusters is designated
        
        for(i in 1:n.group)
        {
            subset <- rownames(cell)[cell$group.for.cluster %in% groups[i]]
            data <- t(scalecounts[,colnames(scalecounts) %in% subset])
            
            if(method == "kernlab"){ # use the kernlab method
                if(nystrom.red == FALSE & is.null(nystrom.sample)){
                    nystrom.sample <- ncol(scalecounts)/6
                }
                sc <- kernlab::specc(data,  # rows are cells
                                     centers = k,
                                     kernel = kernel,
                                     kpar = kpar,
                                     nystrom.red = nystrom.red,
                                     nystrom.sample = nystrom.sample,
                                     mod.sample = mod.sample,
                                     iterations = iterations)
                sampleNames <- colnames(scalecounts[,colnames(scalecounts) %in% subset])
                cell$initial.cluster[match(sampleNames,rownames(cell))] <- sc@.Data
            } else if (method == "kknn"){
                set.seed(12345)
                sc <- kknn::specClust(data = data,
                                      centers = k, nn = 5)
                sampleNames <- colnames(scalecounts[,colnames(scalecounts) %in% subset])
                cell$initial.cluster[match(sampleNames,rownames(cell))] <- sc$cluster
            }
            # update progress bar
            if(display.progress){
                setTxtProgressBar(pb, i)
            }
        }
        
    }  else if (length(k) == n.group) { # n value2 for numbers of clusters are designated. one k for a patient
        
        if (!is.null(names(k)) ) { # match name and name of groups
            k <- k[match(groups, names(k))]
        }
        
        for(i in 1:n.group)
        {
            subset <- rownames(cell)[cell$group.for.cluster %in% groups[i]]
            data <- t(scalecounts[,colnames(scalecounts) %in% subset])
            
            if(method == "kernlab"){
                
                if(nystrom.red == FALSE & is.null(nystrom.sample)){
                    nystrom.sample <- ncol(scalecounts)/6
                }
                
                sc <- kernlab::specc(data,
                                     centers = k[i],
                                     kernel = kernel,
                                     kpar = kpar,
                                     nystrom.red = nystrom.red,
                                     nystrom.sample = nystrom.sample,
                                     mod.sample = mod.sample,
                                     iterations = iterations)
                
                sampleNames <- colnames(scalecounts[,colnames(scalecounts) %in% subset])
                cell$initial.cluster[match(sampleNames,rownames(cell))] <- sc@.Data
            } else if(method == "kknn"){
                sc <- kknn::specClust(data = data,
                                      centers = k[i], nn = 5)
                sampleNames <- colnames(scalecounts[,colnames(scalecounts) %in% subset])
                cell$initial.cluster[match(sampleNames,rownames(cell))] <- sc$cluster
            }
            
            # update progress bar
            if(display.progress){
                setTxtProgressBar(pb, i)
            }
        }
        
    }
    
    
    if(display.progress){
        close(pb)
    }
    cell$initial.cluster <- paste(cell$group.for.cluster, cell$initial.cluster, sep = ".")
    
    return(cell)
}



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
#' 
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

