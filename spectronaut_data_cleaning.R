#'  Takes in a data table in Spectronaut format and removes all proteins 
#'  containing a specified substring (default "Cont_").
#'  Uses grep and targets the $PG.ProteinAccession column
#' @param input.data.table The spectronaut import table to clean
#' @param substring.to.remove (optional) Contaminant substring to look for. Default is "Cont_" 
#' 
Spectronaut_RemoveContaminantProteins <- function(input.data.table, substring.to.remove = "Cont_") {
  
  # Check to make sure $PG.ProteinAccessions exists
  if (!"PG.ProteinGroups" %in% colnames(input.data.table)) {
    stop("\nDataset does not contain column `PG.ProteinGroups`.\nPlease make sure this dataset is in the proper Spectronaut format.")
  }
  
  # Identify rows containing the contaminant substring
  rows_to_remove <- grep(substring.to.remove, input.data.table$PG.ProteinGroups)
  
  # Remove rows using the above vector of row indices
  input.data.table <- input.data.table[-rows_to_remove, ]
  
  return(input.data.table)
}

#'  Takes in a data table in Spectronaut format and thresholds a
#'  specific ptm based on localization score (default phospho, 0.75 respectively).
#'  The frequency filter, if specified, allows you to retain all phospho hits if X number
#' @param input.data.table The spectronaut import table to clean
#' @param retain.unmodified.peptides (optional) Choose to retain unmodified peptides in data table for MsStatsPTM adjustment
#' @param ptm.strings (optional) Modification string to look for. Default is "Phospho (STY)"
#' @param filter.ptm.score (optional) Minimum localization score to accept a PSM. Set to 0 if you want to accept all hits
#' @param filter.ptm.frequency (optional) # of times a PTM must past localization filter to hold onto all PSMs. Set to 1 if you want to accept all PSMs if one passes threshold
#' 
Spectronaut_FilterPTMs <- function(input.data.table, retain.unmodified.peptides = TRUE, ptm.strings = c("Phospho (STY)"), filter.ptm.score = 0.75, filter.ptm.frequency = 1) {
  
  nonmodified.peptides <- input.data.table %>% 
    filter(!grepl(c("Phospho \\(STY\\)"), EG.PrecursorId))

  unfiltered.modified.peptides <- input.data.table %>%
    filter(grepl(c("Phospho \\(STY\\)"), EG.PrecursorId))

  filtered.modified.peptides <- unfiltered.modified.peptides %>%
    group_by(EG.PrecursorId) %>%
    filter(any(EG.PTMAssayProbability >= filter.ptm.score, na.rm = TRUE))

  # Here is where you would filter by frequency, but I don't care about that right now

  # Reassemble modified / nonmodified peptides and return
  if (retain.unmodified.peptides) {
    return(bind_rows(filtered.modified.peptides, nonmodified.peptides))
  } else {
    return(filtered.modified.peptides)
  }
  
  
}