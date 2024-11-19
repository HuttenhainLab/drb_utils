fitPoly.MultiplePowers <- function(data,
                                   polyColumn = "rank.time", otherTerms = c("SUBJECT"), yColumn = "LogIntensities",
                                   powerRange = 1:5, maxSelectedPower = 3, minSelectedPower = 1,
                                   pValueOverall = 0.01, pValueIncrement = 0.05,
                                   start.times = c(0), useAbsoluteLog2FC = TRUE){
  # for debugging...
  #print (unique(data$Protein))
  
  if (length (otherTerms) != 1 | !"SUBJECT" %in% otherTerms )
    stop("other terms other than SUBJECT, or no other terms, is not implemented yet")
  
  # some checks...am I getting enough data to not fail...
  if (nrow(data) < 3){
    return (list(error = sprintf("Too few data points: %s", nrow(data))))
  }
  if (length(unique(data[[polyColumn]])) < 2){
    return (list(error = sprintf("Too few %s values: %s", polyColumn,  length(unique(data[[polyColumn]])))))
  }
  
  # this next check is not needed if otherTerms are all factors
  # if (!nrow(data) > length(unique(data[[polyColumn]]))){ #nrow(unique(data[,c(polyColumn), with = FALSE])) ){
  #   return (list(error = sprintf("No groups in %s have multiple data points: %s groups", polyColumn, nrow(data))))
  # }
  
  # theoretically otherTerms could be for numeric values, but typically we want to treat it as a factor.
  #
  for (term in otherTerms){
    if(any (c("integer", "numeric") %in% class(data[[term]])))
      stop(sprintf ("Term %s is numeric or integer type, expect something factor or character", term))
  }
  
  
  ########################################################
  # some accessory functions -- including them within the function makes for easier parallelization
  ########################################################
  makeFormula <- function(polyColumn, otherTerms, yColumn, power){
    if (power == 0) polyString = "1"
    else polyString = sprintf("poly(%s, %d)", polyColumn, power)
    formString <- sprintf ( "%s ~ %s", yColumn, paste0(c(polyString, otherTerms), collapse=" + "))
    #print(formString)
    as.formula(formString)
  }
  
  
  max.timeDelta <- function(timeDeltas, polyColumn = "rank.time", absolute = TRUE){
    if(grepl("delta", polyColumn)) stop("delta is not allowed in polyColumn name")
    deltaColumns <- grep("delta", colnames(timeDeltas), value=TRUE)
    deltaMatrix <- as.matrix (timeDeltas[,c(polyColumn, deltaColumns), with = FALSE], rownames = polyColumn)
    
    if (!absolute){
      bestDelta <- max(deltaMatrix, na.rm=TRUE)
      maxIdx <- which(deltaMatrix == bestDelta, arr.ind=TRUE)
    }else{
      bestDelta <- max(abs(deltaMatrix), na.rm=TRUE)
      maxIdx <- which(abs(deltaMatrix) == bestDelta, arr.ind=TRUE)
      bestDelta <- deltaMatrix[maxIdx[1,1], maxIdx[1,2]]
    }
    
    bestTime <- timeDeltas[[polyColumn]][maxIdx[1,"row"]] # 1 selects the first returned by which in case of ties
    bestComparison <- deltaColumns[[maxIdx[1,"col"]]] 
    
    return (list(bestTime = bestTime, bestComparison = bestComparison, bestDelta = bestDelta) )
    
  }
  
  # delta versus one or more start times using fitted model estimated marginal means
  deltas.perTime.emm <- function(dt,lmout, start.times = c(0,1), polyColumn = "rank.time"){
    if(grepl("delta", polyColumn)) stop("delta is not allowed in polyColumn name")
    
    if (length(intersect(start.times, dt[[polyColumn]])) == 0) start.times <- min(dt[[polyColumn]])
    
    # build table of all by all rank.times and subjects that are in the data 
    rank.times <- unique(dt[[polyColumn]])
    subjects <- unique(lmout$model$SUBJECT)#unique(dt$SUBJECT)
    newData <- data.table::setnames(data.table::data.table(rank.times), polyColumn)[]
    newData <- newData[, .(SUBJECT = subjects), by = c(polyColumn)]
    
    # calculate expected values
    newData[,prediction := predict (lmout, newData)]
    View(data)
    View(newData)
    # collate actual with imputed for the missing
    newData <- merge(newData, data[, c(polyColumn, otherTerms,  yColumn ), with = FALSE], by = c(polyColumn, otherTerms), all.x=TRUE)
    newData[, actual := newData[[yColumn]] ]
    newData[is.na(actual), actual := prediction]
    
    # and means across subjects
    timeMeans <- newData[,.(meanPrediction = mean(prediction), meanActual = mean(actual)), by = c(polyColumn)]
    #then deltas
    timeMeans[,paste0("delta.", start.times) := lapply (start.times,
                                                        function(t){
                                                          ifelse(t <= timeMeans[[polyColumn]] ,meanActual, NA) - # only compare backwards in time
                                                            timeMeans[timeMeans[[polyColumn]] == t,meanActual]
                                                        }
    )]
    
    return (timeMeans)
  }
  
  
  # get the summary f statistic for the whole model.
  # this includes effects of SUBJECT + 
  # not so useful if there is a strong SUBJECT effect...
  pFSummaryLm <- function (model){
    f <- summary(model)$fstatistic
    pF <- pf(f["value"], f["numdf"], f["dendf"], lower.tail=FALSE)
    return (pF)
  }
  ########################################################
  # /end accessory functions
  ########################################################33
  
  #include 0 for a background model
  powerRange <- union(0, powerRange)
  
  #powers must be integers and sorted
  powerRange <- sort(unique(as.integer(powerRange)))
  # powers must be less than the number of time points
  powerRange <- powerRange[powerRange < length(unique(data[[polyColumn]]))]
  
  # do lms, one for each power in powerRange
  lms <-tryCatch({
    lapply (powerRange, function(power) lm(makeFormula(polyColumn, otherTerms, yColumn, power), data = data))
  }, error = function(err){
    # halt execution, but return the error message.
    print (err)
    return(err)
  })
  
  if("error" %in% class(lms)){
    return (list(error = as.character(lms)))
  }
  
  
  # gather some stats
  # general fit statistics
  pF.summaries <- sapply (lms, pFSummaryLm)
  r.squared <- sapply (lms, function(model)summary(model)$r.squared)
  
  # does incrementing a power help a significant amount?
  if (length(powerRange) > 1){
    increment.anovas <- lapply (2:length(lms), function(i)anova(lms[[i-1]], lms[[i]]))
    pF.increment <- sapply (increment.anovas, function(an)an$`Pr(>F)`[2])
  } else pF.increment <- numeric(0)
  pF.increment <- c(pF.summaries[1], pF.increment)
  
  # does any power (compared to power = 0) give a significant result
  if (length(powerRange) > 1){
    polyColumn.anovas <- lapply (2:length(lms), function(i)anova(lms[[1]], lms[[i]]))
    pF.polyColumn <- sapply (polyColumn.anovas, function(an)an$`Pr(>F)`[2])
  } else pF.polyColumn <- numeric(0)
  pF.polyColumn <- c(pF.summaries[1], pF.polyColumn)
  
  # names(pF.polyColumn) <-  paste("power", powerRange, sep = ".")
  # names(pF.increment) <-  paste("power", powerRange, sep = ".")
  # names(pF.summaries) <- paste("power", powerRange, sep = ".")
  
  # get stats in a convenient named vector
  stats <- data.table::data.table(power = powerRange, pF.summaries, pF.polyColumn, pF.increment, r.squared)
  stats.melt <- data.table::melt(stats, id.vars = "power")
  statsVector <- stats.melt[,value]
  names(statsVector) <- paste(stats.melt$variable, stats.melt$power, sep = ".")
  
  # find the best power.  For now it is the highest power significant in both overall and incremental improvement  
  # bestPower <- suppressWarnings ( stats[power <=  maxSelectedPower &
  #                                         power >= minSelectedPower &
  #                                         pF.polyColumn < pValueOverall &
  #                                         pF.increment < pValueIncrement,
  #                                       max(power)] ) # when the table is empty, will get -Inf and a warning
  # if (bestPower < minSelectedPower) # probably -Inf
  #   bestPower <- minSelectedPower
  
  # or just the best pF.polyColumn 
  data.table::setorder(stats, pF.polyColumn, na.last = TRUE)
  bestPower <- stats[power <=  maxSelectedPower & power >= minSelectedPower,][1,power]
  
  #if no power in desired range, just take bestPower from top of table:
  if (is.na(bestPower)){
    bestPower <- stats[1, power]
  }
  
  statsVector["bestPower"] = bestPower
  statsVector["best.pF.polyColumn"] = stats[power == bestPower, pF.polyColumn]
  
  # get predicted deltas from best lm
  bestLM <- lms[[match(bestPower, powerRange)]]
  dpt <- deltas.perTime.emm(data, bestLM, start.times, polyColumn)
  maxDeltaList <- max.timeDelta(dpt, polyColumn, absolute = useAbsoluteLog2FC)
  meanPredictions <- dpt$meanPrediction
  names(meanPredictions) <- paste0 ("prediction.", dpt[[polyColumn]])
  meanActuals <- dpt$meanActual
  names(meanActuals) <- paste0 ("actualWImputed.", dpt[[polyColumn]])
  
  
  # redo the fit using the best power and without poly, to get meaningful coefficients
  # build dummy variables for power terms...
  for(p in 1:bestPower){
    data[,paste("power", p, sep=".") := .(data[[polyColumn]]^p)]
  }
  bestLM <- lm (as.formula(sprintf("%s ~ %s + %s", 
                                   yColumn, 
                                   paste("power", 1:bestPower, sep=".", collapse="+"),
                                   paste0(otherTerms, collapse=" + "))),
                data = data)
  
  
  statsVector <- c(statsVector, coef(bestLM))
  
  return (c(as.list(statsVector), maxDeltaList, as.list(meanPredictions), as.list(meanActuals)))
  
}

nicePolyFits.fullTable <- function (protQuant, splitColumn = "Protein",
                                    polyColumn = "rankTime", otherTerms = "SUBJECT",
                                    yColumn = "LogIntensities", powerRange = 1:3,
                                    minSelectedPower = 3, ...){
  proteinTables <- split(protQuant, protQuant[[splitColumn]])
  
  
  pf.out <-  pbapply::pblapply(proteinTables, fitPoly.MultiplePowers,
                               polyColumn = polyColumn, otherTerms = otherTerms,
                               yColumn = yColumn, powerRange = powerRange, minSelectedPower = minSelectedPower, ...)
  dt.out <- rbindlist(pf.out,idcol = "Protein", fill = TRUE, use.names = TRUE)
  nicePolyFitOutput(dt.out)
}

