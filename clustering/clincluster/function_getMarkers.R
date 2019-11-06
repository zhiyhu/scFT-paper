getMarkers <- function(matrix, group){
    keep <- rowSums(matrix > 1) > 5
    dge <- edgeR::DGEList(counts = matrix[keep,]) # make a edgeR object
    
    info <- rep("control", ncol(matrix))
    info[group == T] <- "group"
    design <- model.matrix(~ 0 + info)
    v <- voom(dge, design, plot = F)
    fit <- lmFit(v, design) # Linear Model for Series of Arrays
    cont.matrix <- makeContrasts(contrasts = "infogroup-infocontrol",levels=design)
    fit <- contrasts.fit(fit, cont.matrix ) # Linear Model for Series of Arrays
    fit <- eBayes(fit)
    
    marker <- topTable(fit, p.value = 0.05, number = Inf, coef = 1, lfc = 0.6, sort.by = "logFC")
    
    logc <- log1p(calculateCPM(matrix))
    v_expr <- logc[match(rownames(marker), rownames(logc)), info == "group"]
    marker$ratio1 <- rowSums(v_expr > 1)/ncol(v_expr)
    v_expr <- logc[match(rownames(marker), rownames(logc)), info != "group"]
    marker$ratio2 <- rowSums(v_expr > 1)/ncol(v_expr)
    marker$gene <- rownames(marker) 
    
    marker  <-  marker[order(marker$logFC, decreasing = T), ]
    return(marker)
    
}