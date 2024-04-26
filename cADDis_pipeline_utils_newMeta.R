#'  Reformats `metadata` `Well` column from manually typed format to cADDis format
#'  Essentially adds an extra 0 to the beginning of single-digit well numbers if necessary
#'  Example input: `A1`
#'  Example output: `A01`
#' @param input.data.table Freshly imported `metadata` table for cADDis pipeline
Reformat.Metadata.Well.Position <- function(input.data.table) {
  
  # Define a function to reformat each element
  reformat_string <- function(string) {
    
    letter <- gsub("[0-9]", "", string)
    numbers <- gsub("[A-Z]", "", string)
    
    # If there's only one digit, add a leading zero
    if (nchar(numbers) == 1) {
      numbers <- paste0("0", numbers)
    }
    
    # Combine the letter and possibly modified number
    formatted_string <- paste0(letter, numbers)
    
    return(formatted_string)
  }
  
  # Apply the reformatting function to the "Well" column
  input.data.table$Well <- sapply(input.data.table$Well, reformat_string)
  
  return(input.data.table)
}

#'  Adds a `Replicate` column to metadata table
#'  For all unique combinations of `Drug` and `[Conc]`, insert a replicate number 
#'  @param input.data.table `metadata` table for the cADDis pipeline before or 
#'  @param force (optional). Recreates `Replicate` column regardless if it already exists 
Reformat.Metadata.Add.Replicates <- function(input.data.table, force = FALSE) {
  
  if (!force) {
    if ("Replicate" %in% names(input.data.table)) {
      stop("Error: Replicate column already exists in metadata table.\n\nEither remove this function, or set force = TRUE to update Replicate column.")
    }
  }
  
  return.table <- input.data.table %>%
    group_by(Drug, `[Conc]`) %>%
    mutate(Replicate = row_number())
  
  return(return.table)
}

#'  Adds an extra 0 to the beginning of single-digit well numbers if necessary
#'  Example input: `A1`
#'  Example output: `A01`
#' @param input.data.table Freshly imported `metadata` table for cADDis pipeline
Reformat.Metadata.Well.Position <- function(input.data.table) {
  
  # Define a function to reformat each element
  reformat_string <- function(string) {
    
    letter <- gsub("[0-9]", "", string)
    numbers <- gsub("[A-Z]", "", string)
    
    # If there's only one digit, add a leading zero
    if (nchar(numbers) == 1) {
      numbers <- paste0("0", numbers)
    }
    
    # Combine the letter and possibly modified number
    formatted_string <- paste0(letter, numbers)
    
    return(formatted_string)
  }
  
  # Apply the reformatting function to the "Well" column
  input.data.table$Well <- sapply(input.data.table$Well, reformat_string)
  
  return(input.data.table)
}

#'  Takes input cADDis data table and pivots it from wide-to-long format
#'  Replaces RAW*** column headers with cADDis read time points pre-pivot
#'  Also converts time from "#min ##second format to float format
#' @param input.data.table Freshly imported baseline or treatment table for cADDis pipeline
#' @returns Long-format cADDis table containing `Well`, `Measurement.Time`, and `Raw.Reading`
Pivot.cADDIs.Data.Table <- function(input.data.table) {
  
  # Reformatted table to return
  return.table <- input.data.table
  
  return.table <- return.table[,-c(2)]
  
  # Extract time points from the first row, ignoring the first Well" column
  time_points <- return.table[1, -c(1)]
  
  # Remove the first row after extraction(timepoints)
  return.table <- return.table[-c(1), ]
  
  # Rename the column names to the extracted time points
  colnames(return.table)[-c(1)] <- time_points
  
  # Melt the data into long format
  return.table <- return.table %>%
    pivot_longer(cols = -Well, names_to = "Measurement.Time", values_to = "Raw.Reading")
  
  # Convert measurement time to numeric
  return.table$Measurement.Time <- sapply(return.table$Measurement.Time, Format.cADDis.Read.Time)
  
  # Convert fluoresence reading to numeric
  return.table$Raw.Reading <- as.numeric(return.table$Raw.Reading)
    
  return(return.table)
}

#'  Formats cADDis time value from string "# min ## s" to float
#'  e.g. "1 min 30 s" will be converted to 1.5
#'  @param input.data.string String time value in the format "# min" or "# min ## s"
Format.cADDis.Read.Time <- function(input.time.string) {
  
  # split string into minute and second components on space characters
  time.components <- strsplit(input.time.string, ' ')[[1]]
  
  # not sure if hours will ever be a measurement, but throw error to let user know
  if (length(time.components) >= 5) {
    stop(paste0("Error: Issue converting cADDis timepoint \"", input.time.string, "\" to floating-point format. Extra components found.\n\n Pipeline has only been set up to handle strings with the form \"# min ## sec\".Double check the format of measurement times in your cADDis baseline & treatment tables."))
  } else if (length(time.components) == 4) {
    minutes <- as.numeric(time.components[1])
    seconds <- as.numeric(time.components[3]) / 60
    return.number <- minutes + seconds
    
    return(return.number)
  } else if (length(time.components) == 2) {
    return.number <- as.numeric(time.components[1])
    return(return.number)
  } else {
    # Throw error if time.components is hasn't fit any expected format so far.
    stop(paste0("Error: Issue converting cADDis timepoint \"", input.time.string, "\" to floating-point format. Unexpected components found.\n\n Pipeline has only been set up to handle strings with the form \"# min ## sec\".Double check the format of measurement times in your cADDis baseline & treatment tables."))
  }
  
  return(return.number)
}

#'  Extracts no-cADDis background from a baseline or treatment data table
#'  @param input.data.table `baseline` or `treatment` data tables to extract background data from
#'  @param no.cADDis.condition Condition label for the no-cADDis wells
#'  @returns A data table containing `Measurement.Time`, `Reading.Mean`, and `Reading.StDev` from only the no.caddis samples
Extract.cADDis.Background.Signal <- function(input.data.table, no.cADDis.condition, target.column = "Raw.Reading") {
  
  if (!(no.cADDis.condition %in% unique(input.data.table$Drug))) {
    stop(paste0("ERROR: The provided condition string \"", no.cADDis.condition, "\" was not found in the input data table. Double check that this condition is spelled correctly and exists in the input dataset."))
  }
  
  # subset no caddis data from the data table
  no.caddis.data <- input.data.table %>%
    filter(Drug %in% no.cADDis.condition)
  
  target.column <- ensym(target.column)
  
  # Calculate mean and standard deviation for each time point
  return.data <- no.caddis.data %>%
    group_by(Measurement.Time) %>%
    summarize(Reading.Mean = mean(!!target.column),
              Reading.StDev = sd(!!target.column))
  
  return(return.data)
}

#'  Corrects raw assay reads using mean read values from no-cADDis samples
#'  @param input.data.table `baseline` or `treatment` data tables to correct. 
#'  @param cADDis.correction Respective no-cADDis background for either baseline or treatment tables
#'  @param conduct.correction TRUE: Subtract no-cADDis background from `input.data.table`. FALSE: Keep `input.data.table` as-is
#'  @returns A background-corrected data table
Background.Correct.Data.Table <- function(input.data.table, cADDis.correction, conduct.correction = TRUE) {
  
  if (!is.null(cADDis.correction) && conduct.correction == TRUE) {
    
    intermediate.data.table <- input.data.table
    
    intermediate.data.table <- merge(input.data.table, cADDis.correction, by = "Measurement.Time")
    
    intermediate.data.table$Corrected.Reading <- intermediate.data.table$Raw.Reading - intermediate.data.table$Reading.Mean
    
    return.data <- intermediate.data.table[, c("Well", "Drug", "[Conc]", "Fluoro", "Replicate", "Measurement.Time", "Raw.Reading", "Reading.Mean", "Corrected.Reading")]
    return(return.data)
    
  } else if (conduct.correction == FALSE) {
    
    intermediate.data.table <- input.data.table
    
    intermediate.data.table <- merge(input.data.table, cADDis.correction, by = "Measurement.Time")
    
    intermediate.data.table$Corrected.Reading <- intermediate.data.table$Raw.Reading
    
    return.data <- intermediate.data.table[, c("Well", "Drug", "[Conc]", "Fluoro", "Replicate", "Measurement.Time", "Raw.Reading", "Reading.Mean", "Corrected.Reading")]
    return(return.data)
    
  } else {
    stop("Error: cADDis.correction cannot be NULL when conduct.correction = TRUE. Make sure cADDis.correction has been defined, or set conduct.correction = FALSE")
  }
  
}

#'  Calculates ΔF/Fₒ ratio per well from the corrected baseline and treatment data
#'  @param corrected.treatment.table No-cADDis-corrected treatment table. 
#'  @param corrected.baseline.table No-cADDis-corrected baseline table. 
#'  @returns A concatenated baseline & treatment data table with ΔF/Fₒ values calculated for all wells and time points
Calculate.Fluoresence.Ratios <- function(corrected.treatment.table, corrected.baseline.table) {
  
  # To calculate ΔF/Fₒ, we need to take the treatment readings, subtract the average baseline reading, and divide by the average baseline reading.
  
  corrected.treatment.table <- results$treatment$corrected.table
  corrected.baseline.table <- results$baseline$corrected.table
  
  # Calculate average baseline reading per well over the entire baseline table
  baseline.summary <- corrected.baseline.table %>%
    group_by(Drug, `[Conc]`, Well) %>%
    summarize(Baseline.Mean = mean(Corrected.Reading),
              Baseline.StDev = sd(Corrected.Reading))
  
  # merge baseline.summary onto the corrected.baseline.table to do the ΔF/Fₒ calculation per well
  merged.baseline.table <- merge(corrected.baseline.table, baseline.summary, by = c("Well"))
  merged.baseline.table$ΔF.Fₒ <- (merged.baseline.table$Corrected.Reading - merged.baseline.table$Baseline.Mean) / merged.baseline.table$Baseline.Mean
  merged.baseline.table$Type <- "Baseline"
  
  # merge baseline.summary onto the corrected.treatment.table to do the ΔF/Fₒ calculation per well
  merged.treatment.table <- merge(corrected.treatment.table, baseline.summary, by = c("Well"))
  merged.treatment.table$ΔF.Fₒ <- (merged.treatment.table$Corrected.Reading - merged.treatment.table$Baseline.Mean) / merged.treatment.table$Baseline.Mean
  merged.treatment.table$Type = "Treatment"
  
  # shift treatment table timepoints by the max time in the baseline table to make a linear plot
  max.baseline.time <- max(merged.baseline.table$Measurement.Time)
  merged.treatment.table$Measurement.Time <- merged.treatment.table$Measurement.Time + max.baseline.time
  
  return.table <- as.data.table(rbind(merged.baseline.table, merged.treatment.table))
  
  # Remove redundant columns from return table
  return.table[, c("Drug.y", "[Conc].y") := NULL]
  
  # rename drug & conc columns
  setnames(return.table, c("Drug.x", "[Conc].x"), c("Drug", "[Conc]"))
  
  return(return.table)
}
