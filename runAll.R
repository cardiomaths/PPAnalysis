# PPAnalysis 

  #Load all source code
  source("R/calcAxisLimits.R")
  source("R/processData.R")
  source("R/filenameSplit.R")
  source("R/fitEC50.R")
  source("R/dirtyCheck.R")

  source("R/processSummary.R")
  source("R/processBarcode.R")
 
  library(plotrix)
  library(ggplot2)
  library(RColorBrewer)
  library(minpack.lm)
  library(tidyverse)
  library(corrplot)
  library(PerformanceAnalytics)
  library(psych)
  
### DS1
  print("Data step 1 processing - extract summary statistics")

  #creates file data/Limits.txt
  calcAxisLimits('data/')  

  processData('data','outputs/CurveFitting','Cohort1','data/singlePointFile.txt')

### DS2
  print("Data step 2 processing - summary reports")
  
  processSummary('outputs/CurveFitting/summaryFileCohort1V1.csv','outputs/SummaryStats', 'Cohort1','Yes')

### DS3
  print("Data step 3 processing - barcodes")
  
  processBarcode('outputs/SummaryStats/summaryFileCohort1V2.csv','outputs/Barcode','Cohort1','Yes')
