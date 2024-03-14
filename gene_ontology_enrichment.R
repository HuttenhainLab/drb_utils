#' @prerequisite clusterProfiler package

#'  Imports current Gene Ontology mappings from BioConductor.
#' @param database.name (optional) GO term source. Default "org.Hs.eg.db".
#' @param ontology.type (optional) Default ontology type. Default "BP".
#' @param key.type (optional) Primary key type. Default is "Uniprot".
#' 
Load.GO.Terms.From.Bioconductor <- function (database.name = "org.Hs.eg.db", ontology.type = "BP", key.type = "UNIPROT"){ 
  
  message ("Using package ", database.name, " version ", packageVersion(database.name))
  
  GO <- clusterProfiler:::get_GO_data(database.name, ontology.type, key.type)

  gmt <- rbindlist(lapply (GO$EXTID2PATHID, function(x) data.table(ont.id = x)), idcol="gene")
  gmt$ont <- GO$PATHID2NAME[gmt$ont.id]
  setcolorder(gmt, c("ont", "gene", "ont.id"))

  return(gmt)
}

#'  Counts the # of terminal branches downstream of a specific
#'  point in a ComplexHeatmap dendrogram.
#' @param dendrogram.branch Starting point to begin counting
#' @returns An integer representing the number of children under the dendrogram branch
#' 
Dendrogram.Cluster.NumRows <- function(dendrogram.branch) {

  if (!is.list(dendrogram.branch)) {
    return(1)  # If branch is a terminal node
  }
  
  return(Dendrogram.Cluster.NumRows(dendrogram.branch[[1]]) + Dendrogram.Cluster.NumRows(dendrogram.branch[[2]]))
}

#'  
#'  Gets the original data matrix indices of the children in the provided
#'  dendrogram cluster.
#' @param dendrogram.branch Cluster to get indices from
#' @returns An numerical vector holding the original data matrix indices of the rows contained in this cluster
#' 
Dendrogram.Cluster.Get.Original.Indices <- function(dendrogram.branch) {

  if (!is.list(dendrogram.branch)) {
    return(dendrogram.branch)  # If branch is a terminal node, return its index
  }
  
  return(c(Dendrogram.Cluster.Get.Original.Indices(dendrogram.branch[[1]]), Dendrogram.Cluster.Get.Original.Indices(dendrogram.branch[[2]])))
}

#'  
#' Takes a table containing GO enrichment results from clusterProfiler::enricher and creates a heatmap
#' of the topN hits per group
#' @param cluster.enrichment.table A long format table containing the clusterProfiler::enricher output for one to multiple clusters. 
#'                                   --> These clusters can be derived from whatever source
#' @param group.column The name of the column containing the group-distinguishing values. Default "Cluster"
#' @param p.value.column  The name of the column containing the p-values used to plot. Default "p.adjust"
#' @param p.value.threshold The maximum p-value allowed before exclusion from downstream analysis
#' @param minimum.term.occurance The minimum # of times a GO term must appear to be included in downstream analysis
#' @param top.n The number of hits per cluster reported. Default = 10
#' 
GO.Enrichment.Multigroup <- function(cluster.enrichment.table, group.column = "Cluster", p.value.column = "p.adjust", p.value.threshold = 0.01, minimum.term.occurance = 1, top.n = 10, ...) {
  
  #check to make sure all variables are defined.
  if (is.null(cluster.enrichment.table) || nrow(cluster.enrichment.table) == 0) {
    stop("Exception: Mandatory argument `cluster.enrichment.table` is either null or contains no rows. Please make sure this variable is populated.")
  }
  if (is.null(group.column) || group.column == "") {
    stop("Exception: Mandatory argument `group.column` is either null or an emtpy string. Please make sure this variable is populated.")
  }
  if (!(group.column %in% colnames(cluster.enrichment.table))) {
    stop("Exception: Argument `group.column` does not match to any column names in `cluster.enrichment.table`. Please double check your spelling, capitalization, etc.")
  }
  if (is.null(p.value.column) || p.value.column == "") {
    stop("Exception: Mandatory argument `p.value.column` is either null or an emtpy string. Please make sure this variable is populated.")
  }
  if (!(p.value.column %in% colnames(cluster.enrichment.table))) {
    stop("Exception: Argument `p.value.column` does not match to any column names in `cluster.enrichment.table`. Please double check your spelling, capitalization, etc.")
  }
  if (is.null(top.n) || !is.numeric(top.n) || top.n <= 0) {
    stop("Exception: Mandatory argument `top.n` is either null, less than zero, or is a non-numeric value.")
  }
  
  # plotting formatting variable
  row_names_gp = gpar(fontsize = 10)
  
  # reorder cluster.enrichment.table by p.value.column
  setorderv(cluster.enrichment.table, cols = p.value.column)
  
  # Filter table for each provided cluster by grabbing the top N hits that are under the supplied pvalue and minimum count threshold. 
  # Next, slice out the top N IDs per cluster
  top.terms.per.bait <- cluster.enrichment.table[get(p.value.column) < p.value.threshold & Count >= minimum.term.occurance, 
                                                .(ID=ID[1:top.n]),
                                                by=group.column
                                              ]
  
  # Get the unique subset of terms across all clusters
  top.terms <- unique(top.terms.per.bait$ID)
  
  # Swap from long to wide format for easier plotting
  main.variable.wide <- dcast (cluster.enrichment.table[ID %in% top.terms], as.formula(paste("Description", group.column, sep="~")), value.var=p.value.column)
  
  # transform data.table to matrix and log10 transform p values.
  main.variable.matrix <- -log10(as.matrix(main.variable.wide, rownames = "Description"))
  
  # set na values to zero.
  main.variable.matrix[is.na(main.variable.matrix)] <- 0
  
  # Special code to fix incorrectly-formatted GO terms in the ID column
  if (all(grepl("^GO_", rownames(main.variable.matrix)))){
    stop("Error: Need to fix MSIGDB GO Names like BP did. Dain, write the code.")
    #rownames(main.mat) <- fixMsigdbGONames(rownames(main.mat))
  }
  
  # Create an identical matrix in parallel, but instead of p values, track the number of times these GO terms appear
  counts.wide <- dcast (cluster.enrichment.table[ID %in% top.terms], as.formula(paste("Description", group.column, sep="~")), value.var="Count")
  counts.matrix <- as.matrix(counts.wide, rownames="Description")
  
  # Create a separate data table which tracks the UniProt gene accessions per cluster
  gene.table <- cluster.enrichment.table[ID %in% top.terms, .(gene = unlist(strsplit(geneID, split="/"))), by = ID]
  
  # optionally clean the names. I don't think this is necessary for now
  # geneTable[,cleanName := fixMsigdbGONames(ID)]
  
  # Takes the GeneRatio column (e.g. '18/53') and removes the numerator and slash from the fraction and saves it to a new column named geneCount
  genes.in.universe.counts <- unique(cluster.enrichment.table[, .( geneCount = as.integer(gsub("[0-9]+/", "", GeneRatio))), by = c(group.column)])
  
  
  if (nrow(genes.in.universe.counts) != length(unique(genes.in.universe.counts[[group.column]]))){
    stop("non-unique gene counts per group. If you didn't combine multiple differently grouped enrichments, this is unexpected. If it is, set annotatePossibleMatches = FALSE")
  }
  columns <- colnames(main.variable.matrix)
  setkeyv(genes.in.universe.counts, group.column)
  
  bottom.barchart <- HeatmapAnnotation(`Group Sizes` = anno_barplot ( genes.in.universe.counts[columns, geneCount] ))
  
  hm <- heatmapNumbered (main.variable.matrix, counts.matrix, NULL, title, max_pAdjust = p.value.threshold, bottom_annotation = bottom.barchart, row_names_gp = row_names_gp,
                         upperThreshold = NULL,...)
  
  invisible(list(geneTable = gene.table, main.mat = main.variable.matrix, counts.mat = counts.matrix, hmList = hm))
}


