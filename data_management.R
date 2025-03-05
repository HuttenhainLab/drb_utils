#'  Generates subdirectories relative to notebook to store intermediate data
#'  and figures generated throughout the pipeline
#' @param data.name A prepended the directory to write out the file to. Default is current working directory
#' 
Create.Pipeline.Directories <- function(data.name = "pipeline") {
  
  data.dir <- paste(data.name, "data", sep="_")
  figure.dir <- paste(data.name, "figures", sep="_")
  
  if (!dir.exists(data.dir)) {
    dir.create(data.dir)
  }
  
  if (!dir.exists(figure.dir)) {
    dir.create(figure.dir)
  }
}

#'  writes out a file prepended with today's date (YYYYMMDD format) as a .csv
#'  for the matching file name with the latest (most recent) matching name
#' @param export.data the data table or frame to export
#' @param filename Filename of the generated .csv
#' @param base.dir (optional) the directory to write out the file to. Default is current working directory
#' 
Save.Csv.With.Timestamp <- function(export.data, filename, base.dir = getwd()) {
  
  # Format the current date as YYYYMMDD
  formatted.date <- format(Sys.Date(), "%Y%m%d")
  date.filename <- paste(formatted.date, filename, sep = "_")
  fwrite(export.data, paste(base.dir, date.filename, sep = "/"))
}


#'  Checks incoming data frame too see if it's in a supported format
#'  for making a contrast matrix for Group Processing
#' @param input.data.frame The MsStats dataframe to check
#' @param type The type of MsStats Dataframe to validate. 
#'             Currently either "" or "PTM".
#'
Check.MsStats.DataProcess.Output <- function (input.data.frame, type = ""){
  if (type == "") {
    if (!"ProteinLevelData" %in% names(input.data.frame) || !"FeatureLevelData" %in% names(input.data.frame) ){
      stop("Error: Input dataframe list does not include ProteinLevelData and FeatureLevelData items.\n\tThis variable should come directly from the `dataProcess` function.")
    }
  } else if (type == "PTM") {
    if (!"PTM" %in% names(input.data.frame) || !"PROTEIN" %in% names(input.data.frame) ){
      stop("Error: Input dataframe list does not include PTM and PROTEIN items.\n\tThis variable should come directly from the `dataSummarizationPTM` function.\n\tIf this is not a PTM experiment, you should be using the standard MsStats workflow.")
    }
  } else {
    stop("DEV ERROR: Wrong type qualifier for wrapper function. Users should not see this error")
  }
}

#'  Generates an MsStats contrast matrix for
#'  Abundance Proteomics MsStats Group Processing
#' @param input.data.frame The MsStats dataframe to check. Should be in the same shape as the data process function output.
#' @param condition.vector A character vector containing the condition-condition comparisons in the format "Condition1-Condition2"
#' @param exact.match Determines if regular expressions or literal direct matches should be used for contrasts
Make.LFQ.Contrast.Matrix <- function(input.data.frame, condition.vector, exact.match = TRUE){
  Check.MsStats.DataProcess.Output(input.data.frame)
  columnNames <- as.character(levels(input.data.frame$ProteinLevelData$GROUP))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if (i != j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastNames <- character(0)
  for (re in condition.vector){
    if (exact.match) {
      contrastNames <- union(contrastNames, grep(paste0("^", re, "$"), contrasts$name, value = TRUE))
    } else {
      contrastNames <- union(contrastNames, grep(re, contrasts$name, value = TRUE))
    }
    
    
  }
  print (contrastNames)
  contrasts <- contrasts[name %in% contrastNames]
  print (contrasts)
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}

#'  Generates an MsStats contrast matrix for 
#'  Post-translational Modification MsStats Group Processing
#' @param input.data.frame The MsStats dataframe to check. Should be in the same shape as the data process function output.
#' @param condition.vector A character vector containing the condition-condition comparisons in the format "Condition1-Condition2"
#' @param exact.match Determines if regular expressions or literal direct matches should be used for contrasts
Check.MsStats.DataProcess.Output <- function(input.data.frame, data.type){
  if (!is.list(input.data.frame)) stop("Input data should be a list.")
  if (!data.type %in% names(input.data.frame)) stop("Specified data type is missing in input data.")
  
  if (data.type == "PTM") {
    if (!("FeatureLevelData" %in% names(input.data.frame$PTM)) || 
        !("ProteinLevelData" %in% names(input.data.frame$PTM))) {
      stop("Error: Input dataframe's PTM list does not include necessary items.")
    }
  }
  # Additional checks for other data types can be added here
}


Make.PTM.Contrast.Matrix <- function(input.data.frame, condition.vector, exact.match = TRUE){
  columnNames <- as.character(levels(input.data.frame$PTM$ProteinLevelData$GROUP))
  
  positives = c()
  negatives = c()
  for (i in columnNames){
    for (j in columnNames){
      if (i != j){
        positives <- c(positives,i)
        negatives <- c(negatives,j)
      }
    }
  }
  contrasts <- data.table (positives, negatives)
  contrasts[,name := paste(positives,negatives, sep="-")]
  
  contrastNames <- character(0)
  for (re in condition.vector){
    if (exact.match) {
      contrastNames <- union(contrastNames, grep(paste0("^", re, "$"), contrasts$name, value = TRUE))
    } else {
      contrastNames <- union(contrastNames, grep(re, contrasts$name, value = TRUE))
    }
    
    
  }
  print (contrastNames)
  contrasts <- contrasts[name %in% contrastNames]
  print (contrasts)
  
  contrastMat <- matrix(rep(0,length(columnNames) * nrow(contrasts)), nrow=nrow(contrasts), ncol=length(columnNames),
                        dimnames = list(contrasts$name, columnNames))
  
  for (i in seq_len(nrow(contrasts))){
    #rows are named according to the positive condition:
    contrastMat[contrasts$name[i],contrasts$negatives[i]] = -1
    contrastMat[contrasts$name[i],contrasts$positives[i]] =  1
  }
  return (contrastMat)  
}
