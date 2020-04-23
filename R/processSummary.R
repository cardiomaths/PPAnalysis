#' Process summary data (step D2),
#'
#' @param inputFile, the pathname of the summary file to be processed
#'
#' @param outputLocation, pathname of directory where generated figures are stored
#'
#' @param  dsname, a cohort name added to summary file name
#'
#' @param singlePointInfo, Yes or No (default), generate data describing response to a single concentration
#'        of an agonist (this requires the data to have been processed in step 1)
#'
#' @return creates figures and a summary file in the outputLocation
#'
#' # processSummary('../outputs/CurveFitting/summaryFileMethodsV1.csv','../outputs/SummaryStats','Methods','Yes')
processSummary <-
  function(inputFile,
           outputLocation,
           dsname = '',
           singlePointInfo = 'No') {
    
    #library(tidyverse)
    #library(corrplot)
    #library(PerformanceAnalytics)
    #library(psych)
    
    df <-
      read.csv(inputFile,
               stringsAsFactors = FALSE,
               header = TRUE)
    
    # Get rid of all dirty data -------
    df <- subset(df, DirtyStatus != "Dirty")
    
    measures <-
      df[!duplicated(df$measure), c("measure")]
    
    outputs <-
      df[!duplicated(df$output), c("output")]
    
    agonists <-
      df[!duplicated(df$agonist), c("agonist")]
    
    df$date <- as.Date(df$date)
    
    # Create summary plots --------
    for (i in 1:length(measures)) {
      #print(measures[i])
      for (p in 1:length(outputs)) {
        #print(outputs[p])
        df.plot <-
          filter(df, measure == measures[i], output == outputs[p])
        
        if (singlePointInfo == "Yes") {
          df.plot.long <-
            gather(df.plot, var, value, Sensitivity, Capacity)
          df.save.long <-
            gather(df.plot,
                   var,
                   value,
                   Sensitivity,
                   Capacity,
                   singlePointResponse)
        } else {
          df.plot.long <- gather(df.plot, var, value, Sensitivity, Capacity)
          df.save.long <-
            gather(df.plot, var, value, Sensitivity, Capacity)
        }
        
        #add quantiles to the dataframe -------
        df.plot.long$upper <-
          ave(
            df.plot.long$value,
            df.plot.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.95, na.rm = TRUE)
          )
        df.plot.long$lower <-
          ave(
            df.plot.long$value,
            df.plot.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.05, na.rm = TRUE)
          )
        df.plot.long$upper2 <-
          ave(
            df.plot.long$value,
            df.plot.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.99, na.rm = TRUE)
          )
        df.plot.long$lower2 <-
          ave(
            df.plot.long$value,
            df.plot.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.01, na.rm = TRUE)
          )
        
        df.plot.cap <- df.plot.long
        df.plot.cap$var <- paste(df.plot.cap$var, "Cap", sep = "")
        df.plot.cap$value <-
          ifelse(df.plot.cap$value > df.plot.cap$upper,
                 df.plot.cap$upper,
                 df.plot.cap$value)
        df.plot.cap$value <-
          ifelse(df.plot.cap$value < df.plot.cap$lower,
                 df.plot.cap$lower,
                 df.plot.cap$value)
        
        df.save.long$upper <-
          ave(
            df.save.long$value,
            df.save.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.95, na.rm = TRUE)
          )
        df.save.long$lower <-
          ave(
            df.save.long$value,
            df.save.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.05, na.rm = TRUE)
          )
        df.save.long$upper2 <-
          ave(
            df.save.long$value,
            df.save.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.99, na.rm = TRUE)
          )
        df.save.long$lower2 <-
          ave(
            df.save.long$value,
            df.save.long[, c("var", "agonist")],
            FUN = function(x)
              quantile(x, 0.01, na.rm = TRUE)
          )
        
        df.save.cap <- df.save.long
        df.save.cap$var <- paste(df.save.cap$var, "Cap", sep = "")
        df.save.cap$value <-
          ifelse(df.save.cap$value > df.save.cap$upper,
                 df.save.cap$upper,
                 df.save.cap$value)
        df.save.cap$value <-
          ifelse(df.save.cap$value < df.save.cap$lower,
                 df.save.cap$lower,
                 df.save.cap$value)
        
        # Create plot to assess assay consistency over time ---------
        plot1 <- ggplot(df.plot.long, aes(x = date, y = value)) +
          geom_point(aes(colour = factor(plate)), size = 0.5) +
          xlab("response") +
          facet_grid(var ~ agonist, scales = "free") +
          ggtitle(paste("Consistency of ", measures[[i]], " for ", outputs[[p]])) +
          theme_bw() +
          theme(
            text = element_text(size = 8),
            axis.text.x = element_text(size = 4),
            strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")),
            legend.box.margin = margin(-13, 0, 0, 0),
            legend.margin = margin(0, 0, 0, 0),
            legend.position = "bottom",
            legend.title = element_blank()
          )
        
        filenam <-
          paste(outputLocation,
                "/Consistency_",
                measures[[i]],
                "_",
                outputs[[p]],
                ".png",
                sep = "")
        
        png(
          filename = filenam,
          width = 15,
          height = 8,
          units = "cm",
          res = 400
        )
        
        print(plot1)
        
        dev.off()

        
        df.tt <- subset(df.plot.long, agonist != "epinephrine") 
        df.tt <- subset(df.tt, agonist != "U46619") 
        # Create plot to assess assay consistency over time ---------
        plot1 <- ggplot(df.tt, aes(x = date, y = value)) +
          geom_point(aes(colour = factor(plate)), size = 0.5) +
          xlab("response") +
          #facet_grid(var ~ agonist) +
          facet_grid(var ~ agonist, scales = "free") +
          ggtitle(paste("Consistency of ", measures[[i]], " for ", outputs[[p]])) +
          theme_bw() +
          theme(
            text = element_text(size = 8),
            axis.text.x = element_text(size = 4),
            strip.text.x = element_text(margin = margin(0.1, 0, 0.1, 0, "cm")),
            legend.box.margin = margin(-13, 0, 0, 0),
            legend.margin = margin(0, 0, 0, 0),
            legend.position = "bottom",
            legend.title = element_blank()
          )
        
        filenam <-
          paste(outputLocation,
                "/ConsistencyPaper_",
                measures[[i]],
                "_",
                outputs[[p]],
                ".png",
                sep = "")
        
        png(
          filename = filenam,
          width = 10,
          height = 8,
          units = "cm",
          res = 400
        )
        
        print(plot1)
        
        dev.off()
        
        # Create Plot demonstrating cohort variation ------------
        plot2 <-
          ggplot(df.plot.long, aes(x = factor(agonist), y = value)) +
          geom_boxplot(outlier.colour = NA) +
          geom_jitter(
            color = "blue",
            alpha = 0.25,
            size = 1,
            width = 0.15
          ) +
          facet_wrap(~ var, scales = "free") +
          geom_text(
            aes(label = ifelse(value > upper2, donor, '')),
            hjust = 1,
            vjust = 0,
            size = 1.5
          ) +
          geom_text(
            aes(label = ifelse(value < lower2, donor, '')),
            hjust = 1,
            vjust = 0,
            size = 1.5
          ) +
          ggtitle(paste("Variation of", outputs[[p]], "(MFI)")) +
          theme_bw() +
          theme(
            plot.title = element_text(
              size = 10,
              family = "Tahoma",
              face = "bold",
              colour = "black"
            ),
            text = element_text(size = 8),
            axis.text.x = element_text(
              size = 6,
              angle = 45,
              hjust = 1
            ),
            strip.text.y = element_text(size = 6, angle = 90),
            legend.position = "none",
            axis.title.x = element_blank()
          )
        
        filenam <-
          paste(outputLocation,
                "/Variation_",
                measures[[i]],
                "_",
                outputs[[p]],
                ".png",
                sep = "")
        
        png(
          filename = filenam,
          width = 15,
          height = 10,
          units = "cm",
          res = 400
        )
        
        print(plot2)
        
        dev.off()
        
        # Create Plot demonstrating cohort distributions ----------
        plot3 <- ggplot(df.plot.long, aes(x = value)) +
          #geom_histogram(bins=20) +
          geom_density(adjust = 0.8, na.rm = TRUE, fill="grey") +
          facet_wrap(factor(agonist) ~ var,
                     scales = "free",
                     ncol = 2) +
          ggtitle(paste("Distributions for ", measures[[i]], "(MFI)")) +
          theme_bw() +
          theme(
            plot.title = element_text(
              size = 8,
              family = "Tahoma",
              face = "bold",
              colour = "black",
              hjust = 0.5
            ),
            strip.background = element_rect(fill = NA, colour = NA),
            text = element_text(size = 8),
            axis.text.x = element_text(
              size = 6,
              angle = 45,
              hjust = 1
            ),
            strip.text.y = element_text(size = 6, angle = 90),
            legend.position = "none",
            axis.title.x = element_blank()
          )
        
        filenam <-
          paste(outputLocation,
                "/Distribution_",
                measures[[i]],
                "_",
                outputs[[p]],
                ".png",
                sep = "")
        
        png(
          filename = filenam,
          width = 9,
          height = 17,
          units = "cm",
          res = 400
        )
        
        print(plot3)
        
        dev.off()
        
        # Create Plot demonstrating cohort correlations (format 1) --------
        df.plot.long$summaryLab <-
          paste(df.plot.long$var, "\n", df.plot.long$agonist, sep = "")
        df.plot.corr <-
          df.plot.long[, c("donor", "summaryLab", "value")]
        df.plot.corr.wide <- spread(df.plot.corr, summaryLab, value)
        cordata = cor(
          as.matrix(df.plot.corr.wide[,-c(1)]),
          use = "complete.obs",
          method = c("spearman")
        )
        
        filenam <-
          paste(outputLocation,
                "/Correlations_",
                measures[[i]],
                "_",
                outputs[[p]],
                ".png",
                sep = "")
        
        png(
          filename = filenam,
          width = 25,
          height = 28,
          units = "cm",
          res = 400
        )
        
        # Create Plot demonstrating cohort correlations (format 2, ellipses)
        corrplot.mixed(
          cordata,
          tl.cex = 0.7,
          lower = "ellipse",
          upper = "number",
          title = paste("Correlations between metrics"),
          mar = c(0, 0, 5, 0),
          tl.offset = 1
        )
        
        dev.off()
        
        #delete columns --------
        df.plot.cap <-
          df.plot.cap[, which(names(df.plot.cap) != "lower")]
        df.plot.cap <-
          df.plot.cap[, which(names(df.plot.cap) != "upper")]
        df.plot.cap <-
          df.plot.cap[, which(names(df.plot.cap) != "lower2")]
        df.plot.cap <-
          df.plot.cap[, which(names(df.plot.cap) != "upper2")]
        df.save.cap <-
          df.save.cap[, which(names(df.save.cap) != "lower")]
        df.save.cap <-
          df.save.cap[, which(names(df.save.cap) != "upper")]
        df.save.cap <-
          df.save.cap[, which(names(df.save.cap) != "lower2")]
        df.save.cap <-
          df.save.cap[, which(names(df.save.cap) != "upper2")]
        
        if (i == 1 && p == 1) {
          df.save <- df.save.cap
        } else
        {
          df.save <- rbind(df.save, df.save.cap)
        }
        
      }
    }
    
    # Save summary file V2 --------------
    df.save.wide <- spread(df.save, var, value)
    df.new <- merge(df, df.save.wide)
    write.csv(df.new,
              file = paste(outputLocation, "/summaryFile", dsname, "V2.csv", sep = ""))
    
    print("finished")
  }