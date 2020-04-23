#' process raw data, Step D1: creates plots for each donor and a summary file
#'
#' @param inputLoc, the full pathname of the data directory
#' 
#' @param outputLoc, pathname to store the results (figures and summary file)
#' 
#' @param dsname, a cohort name added to summary file name
#' 
#' @param singlePointInfo, extract the concentration acheived at a single point and add it to the summary file (this requires a file 'singlePointFile.txt' in the data directory)
#' 
#' @return creates figures and summary files in outputLocation
processData <-
  function(inputLoc,
           outputLoc,
           dsname='',
           singlePointInfo='none') {

    #library(plotrix) #Required to break the plot axis when including vehicle data
    #library(minpack.lm) #Used to fit curve to data

    filenameList <-
      list.files(
        path = inputLoc,
        full.names = TRUE,
        recursive = TRUE,
        pattern = ".csv"
      )

    # load metadata file -----
    wellData <- NULL
    try(wellData <-
          read.csv(
            paste(inputLoc, "/WellData.txt", sep = ""),
            stringsAsFactors = FALSE,
            encoding = "UTF-8"
          ),
        silent = TRUE)
    if (is.null(wellData)) {
      stop(paste("no WellData.txt in ", inputLoc, "/WellData.txt", sep = ""))
    }
    
    # If carrying out single point analysis load the relevant metadata ------
    if (singlePointInfo != "none") {
      print("including single point analysis")
      dfSP <- read.csv(singlePointInfo,
                       stringsAsFactors = FALSE,
                       header = TRUE)
    }
    
    # Extract how many rows to skip when opening datafiles ----
    headerRows <-
      wellData[(wellData$agonist == 'header'), c("cellCol")]
    maxRows <-
      wellData[(wellData$agonist == 'header'), c("cellRow")]
    
    agonists <-
      wellData[!duplicated(wellData$agonist), c("agonist")]
    agonists <-
      agonists[(agonists != "vehicle") &
                 (agonists != "negative") &
                 (agonists != "header") & 
                 (agonists != "beads")]
    
    antibodies <-
      wellData[!duplicated(wellData$antibody), c("antibody")]
    antibodies <-
      antibodies[(antibodies != "") & (antibodies != "isotype")]
    
    inhibitors <-
      wellData[!duplicated(wellData$inhibitor), c("inhibitor")]
    inhibitors <-
      inhibitors[(inhibitors != "") & (inhibitors != "EDTA")]
    
    measures <-
      wellData[!duplicated(wellData$measure), c("measure")]
    measures <-
      measures[(measures != "")]
    
    #Load the ylimits file (used when plotting)
    limitFile <- NULL
    try(limitsFile <-
          read.csv(
            paste(inputLoc, "/limitsFile.txt", sep = ""),
            stringsAsFactors = FALSE,
            encoding = "UTF-8"
          ),
        silent = TRUE)
    if (is.null(limitsFile)) {
      stop("no limitsFile.txt")
    }
    
    #create the structure for the summary file
    summaryFile <-
      matrix(
        nrow = length(filenameList) * length(agonists) * length(measures) * length(inhibitors) * length(antibodies),
        ncol = 20
      )
    summaryFileCount <- 1
    colnames(summaryFile) <-
      c(
        "donor",
        "output",
        "plate",
        "sample",
        "date",
        "inhibitor",
        "measure",
        "agonist",
        "Sensitivity",
        "Hill",
        "fit",
        "Emax",
        "EmaxConcentration",
        "Capacity",
        "minResponse",
        "finalResponse",
        "singlePointResponse",
        "vehicle",
        "DirtyStatus",
        "DirtyError"
      )
    
    # Process the data -----------
    for (i in 1:length(filenameList)) {
      
      # split the filename to extract information --------
      fileNam <-
        filenameSplit(filenameList[i])
      fileDate <- fileNam[1]
      donor <- fileNam[2]
      plate <- fileNam[3]
      control <- fileNam[4]
      sample <- fileNam[5]
      
      if (maxRows > 0) {
        df <-
          read.csv(
            filenameList[i],
            stringsAsFactors = FALSE,
            header = FALSE,
            skip = headerRows,
            nrows = maxRows
          )
      } else {
        df <-
          read.csv(
            filenameList[i],
            stringsAsFactors = FALSE,
            header = FALSE,
            skip = headerRows
          )
      }
      
      print(paste(
        "Donor: ",
        donor,
        ", plate: ",
        plate
      ))
      
      for (k in 1:length(antibodies)) {

        for (l in 1:length(measures)) {

          #Set the label for Y-axis of plots
          if (antibodies[k] == "P-selectin") {
            if (measures[l] == "ppositive") {
              ylabel <- "P-Selectin Exposure (%pos)"
            } else {
              ylabel <- "P-Selectin Exposure (AU)"
            }
          }
          else
          {
            if (measures[l] == "ppositive") {
              ylabel <- "Fibrinogen Binding (%pos)"
            } else {
              ylabel <- "Fibrinogen Binding (AU)"
            }
          }

        for (q in 1:length(inhibitors)) {
          filenam <-
            paste(
              outputLoc,
              '/',
              donor,
              "_",
              plate,
              "_",
              control,
              "_",
              sample,
              "_",
              measures[[l]],
              "_",
              antibodies[[k]],
              "_",
              inhibitors[[q]],
              ".png",
              sep = ""
            )
          
          png(
            filename = filenam,
            width = 4.1,
            height = 5,
            units = "in",
            res = 400
          )
          par(mfrow = c(3, 2))
          par(oma = c(0.2, 0, 0.8, 0.5)) 
          par(mgp = c(1, 0.1, 0)) 
          par(mar = c(2, 2, 1, 1.75)) 
          par(
            ps = 8,
            cex = 1,
            cex.main = 1,
            cex.axis = 0.8,
            cex.lab = 0.8,
            font.main = 1
          )
          
          # extract vehicle wellData info. -----------
          vehicle <-
            wellData[which(
              wellData$agonist == "vehicle" &
                wellData$control != "PPP" &
                wellData$control != "PRP"
            ), c("cellRow", "cellCol", "concentration", "control")]
          
          for (j in 1:length(agonists)) {

            wellRefs <-
              wellData[which((wellData$agonist == agonists[j]) &
                               (wellData$antibody == antibodies[k]) &
                               (wellData$inhibitor == inhibitors[q]) &
                               (wellData$measure == measures[l])
              ), c("cellCol",
                   "cellRow",
                   "concentration",
                   "unit",
                   "control")]
            
            # ddjust the rowno to account for the rows Ive skipped on opening the datafile  ----
            wellRefs$cellRow <- wellRefs$cellRow - headerRows
            
            # extract data, based on metadata ----
            dat <-
              numeric(dim(wellRefs)[1])
            for (m in 1:nrow(wellRefs)) {
              #I'm stripping the data of "%" and ","
              dat[m] <-
                as.numeric(sub("%", "", gsub(",", "", df[wellRefs$cellRow[m], wellRefs$cellCol[m]], fixed =
                                               TRUE)))
            }
            
            # replace NA with 0
            dat[is.na(dat)] <- 0
            
            dfAgonist <- as.data.frame(dat)
            names(dfAgonist) <- "trans"
            agonistConcentrations <- wellRefs["concentration"]
            dfAgonist <- cbind(dfAgonist, agonistConcentrations)
            agonistUnits <- unique(wellRefs["unit"])
            
            # locate the vehicle and extract data -----
              axisControl <-
                unique(wellRefs[, "control"])
              if (axisControl != "none") {
                vehicleRefs <-
                  vehicle[which(vehicle$control == axisControl), c("cellRow", "cellCol")]
                # adjust the rowno to account for the rows Ive skipped on opening the datafile
                vehicleRefs$cellRow <-
                  vehicleRefs$cellRow - headerRows
                # get the vehicle data corresponding to cells defined in WellData.txt
                dfVehicle <-
                  numeric(dim(vehicleRefs)[1])
                for (m in 1:nrow(vehicleRefs)) {
                  dfVehicle[m] <-
                    as.numeric(sub("%", "", gsub(",", "", df[vehicleRefs$cellRow[m], vehicleRefs$cellCol[m]], fixed =
                                                   TRUE)))
                }
                dfVehicle[is.na(dfVehicle)] <- 0
                vehicleMean <- as.numeric(mean(dfVehicle))
                
                left <- min(dfAgonist$concentration) * 0.01
                break1 <-
                  log10(min(dfAgonist$concentration) * 0.1)
                xlabels <-
                  seq(round(min(
                    log10(dfAgonist$concentration)
                  ), digits = 0), round(max(
                    log10(dfAgonist$concentration)
                  ), digits = 0))
                dfAgonist[nrow(dfAgonist) + 1, "concentration"] <-
                  left
                dfAgonist[nrow(dfAgonist), "trans"] <- vehicleMean
              } 
              else {
                xlabels <-
                  seq(round(min(
                    log10(dfAgonist$concentration)
                  ), digits = 0), round(max(
                    log10(dfAgonist$concentration)
                  ), digits = 0))
                vehicleMean <- 0
              }
 
            # transform the x scale so there are no negatives,
            xscaleTrans <- 20
            dfAgonist$concentration <-
              xscaleTrans + log10(dfAgonist$concentration)
            
            # fit the hill curve ----------
            resp <-
              fitEC50(dfAgonist$concentration,
                         dfAgonist$trans,
                         measures[l])
            
            EC50 <- round(resp[1] - xscaleTrans, digits = 3)
            Hill <- resp[2]
            fit <- resp[3]
            
            # check for dirty data ---------
            err <-
              dirtyCheck(
                dfAgonist$concentration,
                dfAgonist$trans,
                EC50,
                fit,
                antibodies[k],
                measures[l],
                vehicleMean,
                dsname
              )
            
            if (singlePointInfo != "none") {
              if (agonists[j] %in% dfSP$agonist &
                  antibodies[k] %in% dfSP$antibody &
                  measures[l] %in% dfSP$measure) {
                dfSP.specificRow <-
                  dfSP[dfSP$agonist == agonists[j] &
                         dfSP$antibody == antibodies[k] &
                         dfSP$measure == measures[l], c("cellCol", "cellRow")]
                dfSP.specificRow$cellRow <-
                  dfSP.specificRow$cellRow - headerRows
                singlePointResponse <-
                  as.numeric(sub("%", "", gsub(",", "", df[dfSP.specificRow$cellRow, dfSP.specificRow$cellCol], fixed =
                                                 TRUE)))
                
                singlePointResonse <- 0
              }
            }
              
              # add a row to the summary file  -----------
              summaryFile[summaryFileCount, 1] <- donor
              summaryFile[summaryFileCount, 2] <- antibodies[k]
              summaryFile[summaryFileCount, 3] <- plate
              summaryFile[summaryFileCount, 4] <- sample
              summaryFile[summaryFileCount, 5] <- fileDate
              summaryFile[summaryFileCount, 6] <- inhibitors[q]
              summaryFile[summaryFileCount, 7] <- measures[l]
              summaryFile[summaryFileCount, 8] <- agonists[j]
              summaryFile[summaryFileCount, 9] <- EC50
              summaryFile[summaryFileCount, 10] <- Hill
              summaryFile[summaryFileCount, 11] <- fit
              summaryFile[summaryFileCount, 12] <-
                max(dfAgonist$trans)
              dfMaxs <-
                dfAgonist[which(dfAgonist$trans > (max(dfAgonist$trans) - max(dfAgonist$trans) *
                                                     0.1)), "concentration"]
              suppressWarnings(EmaxConcentration <-
                round(min(dfMaxs - xscaleTrans), digits = 3))
              summaryFile[summaryFileCount, 13] <- EmaxConcentration
              summaryFile[summaryFileCount, 14] <-
                round(max(dfAgonist$trans) - min(dfAgonist$trans),
                      digits = 2)
              summaryFile[summaryFileCount, 15] <-
                min(dfAgonist$trans)
              summaryFile[summaryFileCount, 16] <-
                round(((
                  dfAgonist$trans[1] - min(dfAgonist$trans)
                ) / (
                  max(dfAgonist$trans) - min(dfAgonist$trans)
                )) * 100, digits = 0)
              if (singlePointInfo != "none") {
                summaryFile[summaryFileCount, 17] <- singlePointResponse
              } else {
                summaryFile[summaryFileCount, 17] <- "none"
              }
              summaryFile[summaryFileCount, 18] <- vehicleMean
              summaryFile[summaryFileCount, 19] <- DirtyStatus
              summaryFile[summaryFileCount, 20] <- Error
              
              summaryFileCount <- summaryFileCount + 1
              
              # plot --------------
              bottom <- min(dfAgonist$trans)
              top <- max(dfAgonist$trans)
              predPoints <-
                seq(
                  from = min(dfAgonist$concentration),
                  to = max(dfAgonist$concentration),
                  by = 0.02
                )
              
              yLow <-
                limitsFile[limitsFile$agonist == agonists[j] &
                             limitsFile$measure == measures[l] &
                             limitsFile$antibody == antibodies[k], "low"]
              yHigh <-
                limitsFile[limitsFile$agonist == agonists[j] &
                             limitsFile$measure == measures[l] &
                             limitsFile$antibody == antibodies[k], "high"]
              
              plot(
                dfAgonist$concentration - xscaleTrans,
                dfAgonist$trans,
                pch = 16,
                col = "red4",
                main = paste(agonists[j]),
                xaxt = "n",
                yaxt = "n",
                ylim = c(yLow, yHigh),
                xlab = paste("log[", agonists[j], "] (", agonistUnits, ")", sep = ""),
                ylab = ylabel
              )
              axis(2, mgp = c(0, .4, 0))
              axis(1,
                   at = xlabels ,
                   labels = xlabels,
                   mgp = c(0, .1, 0))
              if (axisControl != "None") {
                axis.break(1, break1, style = "slash", brw = 0.06)
              }
              
              if (DirtyStatus == "Dirty") {
                xcord <-
                  min(dfAgonist$concentration) + (max(dfAgonist$concentration) - min(dfAgonist$concentration)) /
                  2 - xscaleTrans
                ycord <- yLow + (yHigh - yLow) / 2
                textErr <- "Dirty"
                if (Error == "insufficient response") {
                  textErr <- "Insufficient\nresponse"
                }
                if (Error == "Bad fit") {
                  textErr <- "Dirty\nBad fit"
                }
                if (Error == "EC50 out of allowable range") {
                  textErr <- "Dirty\nOut of range"
                }
                text(
                  xcord,
                  ycord,
                  textErr,
                  cex = 2,
                  col = rgb(0.5, 0.5, 0.5, alpha = 0.4)
                )
              }
              
              if (Hill != 0) {
                lines(
                  predPoints - xscaleTrans,
                  (bottom + (top - bottom) / (1 + (
                    (EC50 + xscaleTrans) / predPoints
                  ) ^ Hill)),
                  col = "black",
                  type = "l"
                )
                
                legend(
                  "topleft",
                  legend = c(paste("EC50 = ", EC50)),
                  bty = "n",
                  cex = 0.8,
                  inset = -0.1
                )
                segments(
                  EC50,
                  yLow,
                  EC50,
                  yHigh,
                  col = "black",
                  lty = 2,
                  xpd = FALSE
                )
              # segments(
              #    EmaxConcentration,
              #    yLow,
              #    EmaxConcentration,
              #    yHigh,
              #    col = "grey",
              #    lty = 2,
              #    xpd = FALSE
              #  )
              }
              
            }
            
            title(
              main = paste(donor, ", ", antibodies[k], ", ", inhibitors[q]),
              outer = T
            )
            
            dev.off()
            
          }
        }
      }
      
    }
    
    summaryFileName <-
      paste(outputLoc, "/summaryFile", dsname, "V1.csv", sep = "")
    
    write.csv(summaryFile, file = summaryFileName, row.names = FALSE)
    
    print("finished")
  }
