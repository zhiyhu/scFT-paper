#' Identify the differentially expressed genes between two populations
#'
#' @description Identify the differentially expressed genes between two populations
#'
#' @param object The \code{SingleCellExperiment} S4 object returned from \code{InitialCluster()}
#' @param idx.1 The index of group 1
#' @param idx.2 The index of group 2
#' @param p.value The p.value to screen the differentially expressed genes and to be passed to the \code{limma::topTable} （default: 0.01)
#' @param low.pct cutoff needed to identify DE genes. The cutoff of percentage of cells expressing the genes.
#' @param log.FC the least log-fold change needed for DE genes.
#' @param print.process print out the progress or not
#'
#' @details Identify the differentially expressed genes between two populations
#'
#' @return The table of differentially expressed genes and results of differential expression analysis
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
#' @seealso limma voom
#'
#' @export
IdentifyDEGenes <- function(object, idx.1, idx.2,
                            p.value = 0.05,
                            low.pct = 0.3, # expr.quantile = 0.25,
                            log.FC = 0.8,
                            print.process = F) {

    matrix <- counts(object)[,c(idx.1, idx.2)] # subset count matrix
    sample.info <- c(rep("group1",length(idx.1)), rep("group2",length(idx.2))) # set sample information
    dge <- edgeR::DGEList(matrix) # make a edgeR object

    design <- model.matrix(~sample.info)
    colnames(design) <- names(table(sample.info))

    if(print.process){
      print("Start DE analysis")
    }

    # limma voom
    v <- limma::voom(dge, design, plot=F)
    fit <- limma::lmFit(v, design) # Fit linear model for each gene given a series of arrays.
    fit <- limma::eBayes(fit)

    # calculate the proportion of cells that express the gene in each group
    index2.1 <- rowSums(matrix[,sample.info == "group1"] > 0)/length(idx.1)
    index2.2 <- rowSums(matrix[,sample.info == "group2"] > 0)/length(idx.2)
    index2 <- rowMax(cbind(index2.1, index2.2)) # the max proportion
    names(index2) <- rownames(matrix)

    # res.all <-  limma::topTable(fit, coef = 2, number = Inf)
    # res.all$ratio <- index2[match(rownames(res.all), names(index2))]
    if(print.process){
      print("Extract DE results")
    }
    res <- limma::topTable(fit, coef = 2, number = Inf, p.value = p.value) # get the DE genes
    res$ratio <- index2[match(rownames(res), names(index2))]
    res$ratio1 <- index2.1[match(rownames(res), names(index2.1))]
    res$ratio2 <- index2.2[match(rownames(res), names(index2.2))]

    res <- res[abs(res$logFC) > log.FC,] # cutoff by logFC
    res <- res[res$ratio > low.pct,] # cutoff by proportion of expression

  return(res)
}

