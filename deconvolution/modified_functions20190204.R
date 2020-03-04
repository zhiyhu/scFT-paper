
## adapted functions
## 2019 Feb 04


bseq_my <- function(markers){ # adapted from the BSEQ-sc, see: https://github.com/shenorrLab/bseqsc
  # Scale the data
  clusters <- "cell.type"
  samples <- "Sample"
  clusters <- xbioc::pVar(fresheset, clusters)
  samples <- xbioc::pVar(fresheset, samples)
  
  # compute total count
  tot <- colSums(exprs(fresheset))
  # compute CPM
  cpm <- sweep(exprs(fresheset), 2L, tot, '/') # divided by the total counts
  # re-scale with within cell type average
  sc <- as.character(paste(clusters, samples))
  tot_map <- sapply(split(tot, sc), mean)
  cpm <- sweep(cpm, 2L, tot_map[sc], '*') # multiple by the average
  
  # result
  res <- fresheset
  exprs(res) <- cpm
  
  # compute averages on markers
  ids <- intersect(unlist(markers), rownames(res))
  x <- cpm[ids, , drop = FALSE]
  clusters <- as.character(pVar(x, clusters))
  # samples <- as.character(pVar(x, samples))
  # within each cell type
  res <- sapply(unique(clusters), function(ct){ # Average of cpm in clusters
    rowMeans(x[, clusters %in% ct, drop = FALSE])
  })
  
  res <- res[, names(markers)]
  res
}

