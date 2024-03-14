Kinase.Activities.Barplot <- function(scores, kinaseMapped, 
                                    max_pValue = 1.0, max_fdr = 1.0, min_N = 2,
                                    sigKinases = NULL, reverse = FALSE,
                                    useMonoFont = FALSE, useViolin = FALSE,
                                    useSEA = FALSE, reorder = TRUE,
                                    ncol = 3, labelPoints = FALSE){
  by.col <- c("CTRL_GENE_NAME")
  
  if ("Label" %in% colnames(kinaseMapped)){
    #expected in kinaseMapped, even with just 1 label
    if ("Label" %in% colnames(scores)){
      # if  we get here, make sure we are joining by Label
      by.col <- c(by.col, "Label")
    } else if(length(unique(kinaseMapped$Label)) != 1){
      # we're going to ignore label, so make sure data is only from one Label
      warning("Multiple Labels detected in kinaseMapped, but no Label information in scores.  You might be combining log2FC from multiple contrasts inappropriately: ", paste0(unique(kinaseMapped$Label), collpase = ", "))
    }
  }
  
  b <- merge (scores, kinaseMapped, by = by.col)
  
  if (is.null(sigKinases)){
    if (useSEA){
      sigKinases <-  unique(scores[pval.sea< max_pValue & padj.sea < max_fdr & size >= min_N]$CTRL_GENE_NAME)
    }else{
      sigKinases <-  unique(scores[pValue< max_pValue & fdr.BH < max_fdr & N >= min_N]$CTRL_GENE_NAME)
    }
  }
  
  if (length(sigKinases) == 0){
    return ("No significant kinases")
  }
  
  if (reorder == TRUE){
    scoreSummary <- scores[,.(meanZ = mean(Z, na.rm = TRUE)), by = CTRL_GENE_NAME]
    setorder(scoreSummary, meanZ)
    if (useSEA){
      scoreSummary <- scores[,.(meanNES = mean(NES, na.rm = TRUE)), by = CTRL_GENE_NAME]
      setorder(scoreSummary, meanNES)
    }
    kinasesSorted <- scoreSummary[,CTRL_GENE_NAME]
  }else kinasesSorted <- sigKinases
  
  b[,CTRL_GENE_NAME := factor(CTRL_GENE_NAME, levels = kinasesSorted)]
  
  if (reverse){
    p <- ggplot (b[CTRL_GENE_NAME %in% sigKinases,], aes(x=log2FC, y = Label, fill = sigScore, col = sigScore, label = phSiteCombo)) + 
      facet_wrap(facets = ~CTRL_GENE_NAME, ncol = ncol)
  }else{
    p <- ggplot (b[CTRL_GENE_NAME %in% sigKinases,], aes(x=log2FC, y = CTRL_GENE_NAME, fill = sigScore, col = sigScore, label = phSiteCombo)) + 
      facet_wrap(facets = ~Label, ncol = ncol)
  }
  if (useSEA){
    p <- p + aes(fill = sigScore.sea, col = sigScore.sea)
  }
  p <- p +
    geom_vline(xintercept=0.0, lty="dotted", col = "black") + 
    geom_jitter(width=0.0, height=0.1, col="black", alpha=0.5)
  if (useViolin){
    p <- p + geom_violin(scale = "area",  alpha=0.7) 
  } else{
    p <-  p + geom_boxplot( varwidth=FALSE, alpha=0.7,outlier.shape = NA)
  }
  p <- p + 
    scale_color_gradient2(low = "blue", high="red", midpoint=0.0) + 
    scale_fill_gradient2(low = "blue", high="red", midpoint=0.0) + 
    theme_bw() 
  
  if (useMonoFont) p <- p + theme(axis.text.y = element_text( size = 10, family = "mono"))
  
  
  if (labelPoints){
    p <- p + ggrepel::geom_text_repel(data = b[CTRL_GENE_NAME %in% sigKinases & abs(sigScore) > 1.5,
                                               .SD[which.max(abs(log2FC))], by = .(Label, CTRL_GENE_NAME)], size = 2, min.segment.length = 0 )
  }
  
  
  
  return (p)
}