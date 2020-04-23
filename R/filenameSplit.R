#' splits filename into its component parts
#'
#' @param filename, 
#' 
#' @return date,donor,plate,control and sample no.
filenameSplit <-
  function(filename) {
    
    filenamsplit1 <- strsplit(filename, "/")
    filenamsplit2 <-
      strsplit(unlist(filenamsplit1)[length(unlist(filenamsplit1))], "_")
    lengthSplit2 <- length(unlist(filenamsplit2))
    
    date <- as.character(as.Date(as.character(unlist(filenamsplit2)[2]), "%y%m%d"))
    donor <- unlist(filenamsplit2)[3]

    if (lengthSplit2 == 3) {
      donor <- unlist(strsplit(unlist(filenamsplit2)[3], "[.]"))[1]
      plate <- "none"
      control <- "none"
      sample <- "1"
    }
    
    if (lengthSplit2 == 4) {
      plate <- unlist(strsplit(unlist(filenamsplit2)[4], "[.]"))[1]
      control <- "none"
      sample <- "1"
    }
    
    if (lengthSplit2 == 5) {
      plate <- unlist(filenamsplit2)[4]
      control <-
        unlist(strsplit(unlist(filenamsplit2)[5], "[.]"))[1]
      sample <- "1"
    }
    
    if (lengthSplit2 == 6) {
      plate <- unlist(filenamsplit2)[4]
      control <- unlist(filenamsplit2)[5]
      sample <-
        unlist(strsplit(unlist(filenamsplit2)[6], "[.]"))[1]
    }
    
    return(c(date,donor,plate,control,sample))

  }
