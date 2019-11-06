#' Initially cluster within groups by kmeans
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
# InitialCluster_spectral <- function(scalecounts, cell, k, display.progress = T,
#                                     kernel = "rbfdot", kpar = "automatic",
#                                     nystrom.red = FALSE, nystrom.sample = NULL,
#                                     mod.sample = 0.85, iterations = 200,
#                                     method = "kernlab",
#                                     ...){
#     # subset <- rownames(cell)[cell$group.for.cluster %in% groups[i]]
#     # data <- t(scalecounts[,colnames(scalecounts) %in% subset])
# 
#     sc <- kknn::specClust(data = t(scalecounts),
#                           centers = k, nn = 5)
#     sampleNames <- colnames(scalecounts[,colnames(scalecounts) %in% subset])
#     cell$initial.cluster[match(sampleNames,rownames(cell))] <- sc$cluster
#   
#     return(cell)
# }

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

  if (length(k) != 1 & length(k) != n.group ) { # wrong length of k

    print("Incorrect length of k. Please check...")
    return(cell)

  } else if (length(k) == 1) { # only a uniform value for number of clusters is designated

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

#' Initially cluster within groups by kmeans
#' @importFrom apcluster apcluster
#' @importFrom apcluster negDistMat
InitialCluster_AP <- function(x, cell){
  groups = unique(cell$group)
  n.group = length(groups) # number of patients
  cell$initial.cluster <- NA
  total <- n.group
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for(i in 1:n.group)
  {
    subset <- rownames(cell)[cell$group %in% groups[i]]
    apres <- apcluster::apcluster(apcluster::negDistMat(r=2), x = as.matrix(t(x[,colnames(x) %in% subset])), details = F)
    cluster.i <- c()
    for(n in 1:length(apres@clusters))
    {
      cluster.i <- c(cluster.i, rep(n, length(apres@clusters[[n]])))
    }
    cell$initial.cluster[match(names(unlist(apres@clusters)), cell$Sample)] <- cluster.i
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cell$initial.cluster <- paste(cell$group, cell$initial.cluster, sep = ".")
  return(cell)
}

#' Initially cluster within groups by kmeans
InitialCluster_kmeans <- function(x, cell, k){
  groups = unique(cell$group)
  n.group = length(groups) # number of patients
  k.cluster <- k
  cell$initial.cluster <- NA

  total <- n.group
  # create progress bar
  pb <- txtProgressBar(min = 0, max = total, style = 3)
  for(i in 1:n.group) # do clustering for each patient
  {
    subset <- rownames(cell)[cell$group %in% groups[i]]
    km.out <- kmeans(t(x[,colnames(x) %in% subset]), k.cluster, nstart = 20)
    cell$initial.cluster[match(names(km.out$cluster),colnames(x))] <- km.out$cluster
    # update progress bar
    setTxtProgressBar(pb, i)
  }
  close(pb)
  cell$initial.cluster <- paste(cell$group, cell$initial.cluster, sep = ".")
  return(cell)
}

#' Differential expression analysis for a pair of clusters
#'
#' @description Run the differential expression analysis between a pair of clusters
#'
#' @param x A vector that defines which pair of clusters should be compared here
#' @param fit An MArrayLM object containing the result of the fits. It is returned from \code{limma:lmFit}
#' @param p.value The p.value to screen the differentially expressed genes and to be passed to the \code{limma::topTable}
#' @param length.contrast The length of contrast for the \code{limma::contrast.fit}
#'
#' @return The number of differentially expressed genes between this pair of clusters.
#'
#' @importFrom limma voom
#' @importFrom limma lmFit
#' @importFrom limma contrasts.fit
#' @importFrom limma eBayes
#' @importFrom limma topTable
#'
#' @references \code{citation("limma")}
#'
RunLimmaVoom <- function(x, fit, p.value, length.contrast)
{
  design <- model.matrix(~initial.cluster)
  colnames(design) <- names(table(initial.cluster))

  v <- limma::voom(dge, design, plot=F)
  fit <- limma::lmFit(v, design) # Fit linear model for each gene given a series of arrays.
  # return MArrayLM fitted model object

  i <- x[1]
  j <- x[2]
  contrast <- rep(0, length.contrast)
  contrast[j] <- 1
  if(i > 1)
  {
    contrast[i] <- -1
  }
  fit2 <- limma::contrasts.fit(fit, contrast)
  fit2 <- limma::eBayes(fit2)
  res <- limma::topTable(fit2, number = Inf, p.value = p.value)

  return(as.numeric(nrow(res)))
}
