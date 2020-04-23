#' calculate the Y-limits for use when plotting the data
#'
#' @param inputLocation, the full pathname of the data directory
#'
#' @param assay, the name you wish assigned to the set of data you are analysis,
#'
#' @return create file limisFile.txt that will be used to set Y limits in figures
calcAxisLimits <-
  function(inputLocation, assay = '') {
    print(paste("calculating axis limits for: ", inputLocation))
    
    filenameList <-
      list.files(
        path = inputLocation,
        full.names = TRUE,
        recursive = TRUE,
        pattern = ".csv"
      )
    
    excludedList <-
      list.files(
        path = inputLocation,
        full.names = TRUE,
        recursive = TRUE,
        pattern = "info|Archive|exclude"
      )
    filenameList <- filenameList[!filenameList %in% excludedList]
    
    # load default wellData file ----
    wellData <- NULL
    try(wellData <-
          read.csv(
            paste(inputLocation, "/WellData.txt", sep = ""),
            stringsAsFactors = FALSE,
            encoding = "UTF-8"
          ),
        silent = TRUE)
    if (is.null(wellData)) {
      stop("no wellData.txt")
    }
    wellDataOrig <- wellData
    
    # extract number of rows to skip when opening datafiles
    headerRows <-
      wellData[(wellData$agonist == 'header'), c("cellCol")]
    
    agonists <-
      wellData[!duplicated(wellData$agonist), c("agonist")]
    agonists <-
      agonists[(agonists != "vehicle") &
                 (agonists != "negative") & (agonists != "header")]
    antibodies <-
      wellData[!duplicated(wellData$antibody), c("antibody")]
    antibodies <-
      antibodies[(antibodies != "") & (antibodies != "isotype")]
    measures <-
      wellData[!duplicated(wellData$measure), c("measure")]
    measures <-
      measures[(measures != "")]
    
    lenAntibodies <- length(antibodies)
    if (lenAntibodies == 'NA') {
      lenAntibodies <- 1
    }
    
    # create the limits file structure
    limitsFile <-
      data.frame(matrix(
        nrow = length(agonists) * length(measures) * lenAntibodies,
        ncol = 5
      ))
    limitsFileCount <- 1
    colnames(limitsFile) <-
      c("agonist",
        "antibody",
        "measure",
        "high",
        "low")
    limitsFile$high <-
      rep(100, length(agonists) * length(measures) * lenAntibodies)
    limitsFile$low <-
      rep(0, length(agonists) * length(measures) * lenAntibodies)
    fileRow <- 1
    
    for (k in 1:length(antibodies)) {
      for (l in 1:length(measures)) {
        for (i in 1:length(agonists)) {
          limitsFile[fileRow, "antibody"] <- antibodies[k]
          limitsFile[fileRow, "agonist"] <- agonists[i]
          limitsFile[fileRow, "measure"] <- measures[l]
          fileRow <- fileRow + 1
        }
      }
    }
    
    for (i in 1:length(filenameList)) {
      filenam <- filenameList[[i]]
      
      df <-
        read.csv(
          filenameList[i],
          stringsAsFactors = FALSE,
          header = FALSE,
          skip = headerRows
        )
      
      for (k in 1:length(antibodies)) {
        for (l in 1:length(measures)) {
          # ppositive has a max of 100 so this doesn't need altering
          if (measures[l] != "ppositive") {
            for (j in 1:length(agonists)) {
              wellRefs <-
                wellData[which((wellData$agonist == agonists[j]) &
                                 (wellData$antibody == antibodies[k]) &
                                 (wellData$measure == measures[l])
                ), c("cellCol", "cellRow", "concentration", "unit")]
              
              
              # adjust the rowno to account for the rows I've skipped on opening the datafile
              wellRefs$cellRow <- wellRefs$cellRow - headerRows
              
              # get the data corresponding to cells defined in WellData.txt
              dfAgonist <-
                numeric(dim(wellRefs)[1])
              
              for (m in 1:nrow(wellRefs)) {
                dfAgonist[m] <-
                  as.numeric(gsub(",", "", df[wellRefs$cellRow[m], wellRefs$cellCol[m]], fixed =
                                    TRUE))
              }
              dfAgonist[is.na(dfAgonist)] <- 0
              
              if (assay == 'PBA') {
                yHigh <- max((10 ^ (-dfAgonist)) * 100)
                yLow <- min((10 ^ (-dfAgonist)) * 100)
              } else {
                yHigh <- max(dfAgonist)
                yLow <- min(dfAgonist)
              }
              
              if (yHigh > limitsFile[limitsFile$agonist == agonists[j] &
                                     limitsFile$measure == measures[l] &
                                     limitsFile$antibody == antibodies[k], "high"]) {
                limitsFile[limitsFile$agonist == agonists[j] &
                             limitsFile$measure == measures[l] &
                             limitsFile$antibody == antibodies[k], "high"] <-
                  yHigh
              }
              if (yLow < limitsFile[limitsFile$agonist == agonists[j] &
                                    limitsFile$measure == measures[l] &
                                    limitsFile$antibody == antibodies[k], "low"]) {
                limitsFile[limitsFile$agonist == agonists[j] &
                             limitsFile$measure == measures[l] &
                             limitsFile$antibody == antibodies[k], "low"] <-
                  yLow
              }
            }
          }
          
        }
        
      }
    }
    
    limitsFileName <-
      paste(inputLocation, "/limitsFile.txt", sep = "")
    
    write.csv(limitsFile, file = limitsFileName, row.names = FALSE)
    
    print("finished")
  }
