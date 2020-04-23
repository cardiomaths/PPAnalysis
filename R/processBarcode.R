#' Function to produce a barcode for each donor (and a summary csv file),
#'
#' @param inputName, the full pathname of the file that stores summary data
#'
#' @param outputLocation, pathname where generated figures are created
#'
#' @param dsname, a cohort name that will be added to summary file names
#'
#' @param singlePointInfo, Yes or No (default), generate data describing response to a single concentration
#'        of an agonist (this requires the data to have been processed in step 1 and 2)
#'
#' @return create figures in outputLocation
processBarcode <-
  function(inputFile,
           outputLocation,
           dsname = '',
           singlePointInfo = 'No')
  {
    #library(ggplot2)
    #library(RColorBrewer)
    
    print(paste("processBarcode: ", outputLocation, "/summaryFile", dsname, "V3.csv", sep = ""))
    
    df <-
      read.csv(inputFile,
               stringsAsFactors = FALSE,
               header = TRUE)
    
    df <- subset(df, DirtyStatus != "Dirty" | is.na(DirtyStatus))
    df <- subset(df, inhibitor == "vehicle")
    
    df <- subset(df, measure == "MedianFI")
    if (singlePointInfo == "No") {
      df <-
        df[, c("donor",
               "output",
               "agonist",
               "SensitivityCap",
               "CapacityCap")]
      names(df) <-
        c("donor", "output", "agonist", "Sensitivity", "Capacity")
    } else {
      df <-
        df[, c(
          "donor",
          "output",
          "agonist",
          "SensitivityCap",
          "CapacityCap",
          "singlePointResponseCap"
        )]
      names(df) <-
        c("donor",
          "output",
          "agonist",
          "Sensitivity",
          "Capacity",
          "singlePoint")
    }
    
    outputs <-
      df[!duplicated(df$output), c("output")]
    
    agonists <-
      df[!duplicated(df$agonist), c("agonist")]
    
    donors <-
      df[!duplicated(df$donor), c("donor")]
    
    #create structure for the summary file
    summaryFile <-
      matrix(nrow = length(donors) * length(agonists) * length(outputs),
             ncol = 17)
    summaryFileCount <- 1
    colnames(summaryFile) <-
      c(
        "donor",
        "output",
        "agonist",
        "Sensitivity",
        "Capacity",
        "kvalue",
        "metric",
        "donorMetric",
        "donorSensitivity",
        "donorCapacity",
        "sensitivitymin",
        "Sensitivitymax",
        "Capacitymin",
        "Capacitymax",
        "singlePoint",
        "singlePointmin",
        "singlePointmax"
      )
    
    df[, "Sensitivitymin"] <- 0
    df[, "Sensitivitymax"] <- 0
    df[, "Capacitymin"] <- 0
    df[, "Capacitymax"] <- 0
    df[, "singlePointmin"] <- 0
    df[, "singlePointmax"] <- 0
    
    # Calculate the min and max of EC50 and maxGain, for use in scalings
    for (i in 1:length(outputs)) {
      for (j in 1:length(agonists)) {
        Sensitivitymin <-
          min(df$Sensitivity[which(df$agonist == agonists[j] &
                                     df$output == outputs[i])], na.rm = TRUE)
        Capacitymin <-
          min(df$Capacity[which(df$agonist == agonists[j] &
                                  df$output == outputs[i])], na.rm = TRUE)
        Sensitivitymax <-
          max(df$Sensitivity[which(df$agonist == agonists[j] &
                                     df$output == outputs[i])], na.rm = TRUE)
        Capacitymax <-
          max(df$Capacity[which(df$agonist == agonists[j] &
                                  df$output == outputs[i])], na.rm = TRUE)
        df$Sensitivitymin[which(df$agonist == agonists[j] &
                                  df$output == outputs[i])] <-
          Sensitivitymin
        df$Capacitymin[which(df$agonist == agonists[j] &
                               df$output == outputs[i])] <-
          Capacitymin
        df$Sensitivitymax[which(df$agonist == agonists[j] &
                                  df$output == outputs[i])] <-
          Sensitivitymax
        df$Capacitymax[which(df$agonist == agonists[j] &
                               df$output == outputs[i])] <-
          Capacitymax
        
        if (singlePointInfo != "No")  {
          singlePointmin <-
            min(df$singlePoint[which(df$agonist == agonists[j] &
                                       df$output == outputs[i])], na.rm = TRUE)
          singlePointmax <-
            max(df$singlePoint[which(df$agonist == agonists[j] &
                                       df$output == outputs[i])], na.rm = TRUE)
          df$singlePointmin[which(df$agonist == agonists[j] &
                                    df$output == outputs[i])] <-
            singlePointmin
          df$singlePointmax[which(df$agonist == agonists[j] &
                                    df$output == outputs[i])] <-
            singlePointmax
        }
      }
    }
    
    df$Sensitivity <-
      (1 - (df$Sensitivity - df$Sensitivitymin) / (df$Sensitivitymax - df$Sensitivitymin))
    
    df$Capacity <-
      (df$Capacity - df$Capacitymin) / (df$Capacitymax - df$Capacitymin)
    
    df$metric <-
      round((df$Capacity + df$Sensitivity) / 2, digits = 2)
    
    if (singlePointInfo != "No")  {
      df$sPoint <-
        (df$singlePoint - df$singlePointmin) / (df$singlePointmax - df$singlePointmin)
    }
    
    #Loop through the donor data - outputs and agonists and calculate the final metric
    for (k in 1:length(donors)) {
      dfdonor <- subset(df, donor == donors[k])
      
      #This calculates a donors aggregate response metric
      tempResponse <- 0
      tempSensitivity <- 0
      tempCapacity <- 0
      tempN <- 0
      for (i in 1:length(outputs)) {
        dfdonor.output <- subset(dfdonor, output == outputs[i])
        
        for (j in 1:length(agonists)) {
          dfdonor.agonist <-
            subset(dfdonor.output, agonist == agonists[j])
          
          if (dim(dfdonor.agonist)[1] != 0) {
            tempResponse <-
              tempResponse + dfdonor.agonist$metric
            tempSensitivity <-
              tempSensitivity + dfdonor.agonist$Sensitivity
            tempCapacity <-
              tempCapacity + dfdonor.agonist$Capacity
            tempN <- tempN + 1
          }
        }
      }
      
      dfdonor[, "donorMetric"] <-
        round(tempResponse / tempN, digits = 2)
      
      donorMetric <- round(tempResponse / tempN, digits = 2)
      
      dfdonor[, "donorSensitivity"] <-
        round(tempSensitivity / tempN, digits = 2)
      
      donorSensitivity <- round(tempSensitivity / tempN, digits = 2)
      
      dfdonor[, "donorCapacity"] <-
        round(tempCapacity / tempN, digits = 2)
      donorCapacity <- round(tempCapacity / tempN, digits = 2)
      
      
      #This is a work around... donors with high sensitivity but low capacity are very wide ... or vis versa,
      # they won't fit within x or y axis limits (1 - 2)
      #so, for plotting only I don't let capacity or sensitivity be too low (near 0) and restrict the maximum
      #sensitivity * kvalue (or capacity*kvalue) to be within the plotting limits (<2)
      dfdonorTemp <- dfdonor
      
      dfdonorTemp <-
        transform(dfdonorTemp,
                  Sensitivity = ifelse(Sensitivity < 0.1, 0.1, Sensitivity))
      dfdonorTemp <-
        transform(dfdonorTemp, Capacity = ifelse(Capacity < 0.1, 0.1, Capacity))
      
      dfdonorTemp$kvalue <-
        sqrt((dfdonorTemp$Sensitivity + dfdonorTemp$Capacity) / (2 * (
          dfdonorTemp$Sensitivity * dfdonorTemp$Capacity
        )))
      
      dfdonorTemp <-
        transform(dfdonorTemp,
                  Sensitivity = ifelse(
                    Sensitivity * kvalue > 2,
                    round(1.9 /
                            kvalue, digits = 2),
                    Sensitivity
                  ))
      dfdonorTemp <-
        transform(dfdonorTemp,
                  Capacity = ifelse(Capacity * kvalue > 2, round(1.9 / kvalue, digits = 2), Capacity))
      
      grid.xmax <- 2
      grid.ymax <- 2
      
      dfplot <- dfdonorTemp
      dfplot <- subset(dfplot, agonist != "epinephrine")
      dfplot <- subset(dfplot, agonist != "U46619")
      
      bb <-
        ggplot(
          dfplot,
          aes(
            xmin = (grid.xmax - kvalue * Sensitivity) / 2,
            xmax = kvalue * Sensitivity + (grid.xmax - kvalue * Sensitivity) / 2,
            ymin = (grid.ymax - kvalue * Capacity) / 2,
            ymax = kvalue * Capacity + (grid.ymax - kvalue * Capacity) / 2
          )
        ) +
        geom_rect(aes(fill = metric)) +
        #ggtitle(paste(donors[k], ", response metric = ", dfdonorTemp$donorMetric)) +
        ggtitle(paste(donors[k])) +
        scale_x_continuous(name = "Sensitivity", limits = c(0, grid.xmax)) +
        scale_y_continuous(name = "Capacity", limits = c(0, grid.ymax)) +
        # geom_text(aes(label = metric),
        #           x = 1,
        #           y = -0.02,
        #           size = 1.5) +
        scale_fill_gradientn(colors = brewer.pal(9, "Greys")[c(3, 4, 5, 6, 7, 8, 9)]) +
        facet_wrap(output ~ agonist, nrow = 1) +
        #theme(plot.background=element_rect(fill="darkred")) +
        #theme_bw() +
        theme(
          plot.background = element_rect(fill = "white"),
          plot.margin = unit(c(0.01, 0.1, 0.01, 0.01), "cm"),
          text = element_text(size = 8),
          plot.title = element_text(
            size = 6,
            colour = "black",
            vjust = -1
          ),
          axis.title.x = element_text(
            colour = "black",
            vjust = 2,
            margin = margin(0.01, 0, 0.01, 0, "cm")
          ),
          axis.title.y = element_text(colour = "black", vjust = -1),
          strip.background = element_rect(fill = "gray96"),
          strip.text = element_text(size = rel(0.5), colour = 'darkred'),
          strip.text.x = element_text(margin = margin(0.03, 0, 0.03, 0, "cm")),
          strip.text.y = element_text(margin = margin(0, 0.03, 0, 0.03, "cm")),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          panel.spacing = unit(0, "lines"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "none"
        )
      
      filenam <-
        paste(outputLocation,
              "/Barcode_",
              donors[k],
              ".png",
              sep = "")
      
      png(
        filename = filenam,
        width = 5,
        height = 2,
        units = "cm",
        res = 400
      )
      
      print(bb)
      
      dev.off()
      
      outputs <-
        dfdonor[!duplicated(dfdonor$output), c("output")]
      
      #Create summary file one row per donor/output/agonist
      for (i in 1:length(outputs)) {
        dfdonor.output <- subset(dfdonor, output == outputs[i])
        dfdonorTemp.output <-
          subset(dfdonorTemp, output == outputs[i])
        
        agonists <-
          dfdonor.output[!duplicated(dfdonor.output$agonist), c("agonist")]
        
        for (j in 1:length(agonists)) {
          dfdonor.agonist <- subset(dfdonor.output, agonist == agonists[j])
          dfdonorTemp.agonist <-
            subset(dfdonorTemp.output, agonist == agonists[j])
          
          #Add to summary file
          summaryFile[summaryFileCount, 1] <- donors[k]
          summaryFile[summaryFileCount, 2] <- outputs[i]
          summaryFile[summaryFileCount, 3] <- agonists[j]
          summaryFile[summaryFileCount, 8] <- donorMetric
          summaryFile[summaryFileCount, 9] <- donorSensitivity
          summaryFile[summaryFileCount, 10] <- donorCapacity
          if (dim(dfdonor.agonist)[1] == 0) {
            summaryFile[summaryFileCount, 4] <- 0
            summaryFile[summaryFileCount, 5] <- 0
            summaryFile[summaryFileCount, 6] <- 0
            summaryFile[summaryFileCount, 7] <- 0
            summaryFile[summaryFileCount, 11] <- 0
            summaryFile[summaryFileCount, 12] <- 0
            summaryFile[summaryFileCount, 13] <- 0
            summaryFile[summaryFileCount, 14] <- 0
            summaryFile[summaryFileCount, 15] <- 0
            summaryFile[summaryFileCount, 16] <- 0
            summaryFile[summaryFileCount, 17] <- 0
          } else {
            summaryFile[summaryFileCount, 4] <-
              dfdonor.agonist$Sensitivity
            summaryFile[summaryFileCount, 5] <-
              dfdonor.agonist$Capacity
            summaryFile[summaryFileCount, 6] <-
              dfdonorTemp.agonist$kvalue
            summaryFile[summaryFileCount, 7] <-
              dfdonor.agonist$metric
            summaryFile[summaryFileCount, 11] <-
              dfdonor.agonist$Sensitivitymin
            summaryFile[summaryFileCount, 12] <-
              dfdonor.agonist$Sensitivitymax
            summaryFile[summaryFileCount, 13] <-
              dfdonor.agonist$Capacitymin
            summaryFile[summaryFileCount, 14] <-
              dfdonor.agonist$Capacitymax
            if (singlePointInfo == "No") {
              summaryFile[summaryFileCount, 15] <- 0
              summaryFile[summaryFileCount, 16] <- 0
              summaryFile[summaryFileCount, 17] <- 0
            } else {
              summaryFile[summaryFileCount, 15] <-
                dfdonor.agonist$sPoint
              summaryFile[summaryFileCount, 16] <-
                dfdonor.agonist$singlePointmin
              summaryFile[summaryFileCount, 17] <-
                dfdonor.agonist$singlePointmax
            }
          }
          summaryFileCount <- summaryFileCount + 1
        }
      }
      
    }
    
    summaryFileName <-
      paste(outputLocation, "/summaryFile", dsname, "V3.csv", sep = "")
    
    write.csv(summaryFile, file = summaryFileName, row.names = FALSE)
    
    print("finished")
    
  }
