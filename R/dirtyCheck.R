#' determine if the data or the fitting process is outside of expected ranges
#' 
#' @param L, data (x-axis)
#' 
#' @param Y, data (y-axis)
#' 
#' @param EC50, 
#' 
#' @param goodness, 
#' 
#' @param antibody, 
#' 
#' @param measure, 
#' 
#' @param vehicle, 
#' 
#' @param dsname
#' 
#' @return 
dirtyCheck <- function(L, Y, EC50, goodness, antibody, measure, vehicle, dsname) {
  DirtyStatus <<- ''
  Error <<- ''

  bottom <- min(Y)
  top <- max(Y)

  rise <- ((top - bottom)/top)*100
  gain <- (top - bottom)
  
  if (EC50 < -9) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'EC50 out of allowable range'
  }
  
  if (EC50 > -3) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'EC50 out of allowable range'
  }
  
  if (goodness == 0) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'Bad fit'
  }

   if ((antibody == 'fibrinogen' || antibody == 'P-selectin') && measure == 'ppositive' && goodness > 400 ) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'Bad fit'
  }

  if (antibody == 'fibrinogen' && measure != 'ppositive' && goodness > 1900 ) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'Bad fit'
  }

  if (antibody == 'P-selectin'  && measure != 'ppositive' && goodness > 4e+7 ) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'Bad fit'
  }
  
  if (vehicle > bottom && (vehicle - bottom)/gain > 0.5) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'Vehicle too high'
  }

if ((antibody == 'fibrinogen' || antibody == 'P-selectin') && measure == 'ppositive' ) {
  if (gain < 10 ) {
  DirtyStatus <<- 'Dirty'
  Error <<- 'insufficient response'
  }
}

try(if (antibody == 'fibrinogen' && measure != 'ppositive' && dsname == 'Methods' &&  gain < 1000) {
  DirtyStatus <<- 'Dirty'
  Error <<- 'insufficient response'
},
silent = TRUE) 
  
try(if (antibody == 'P-selectin' && measure != 'ppositive' && dsname == 'Methods' && gain < 2500) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'insufficient response'
  },
  silent = TRUE) 

  try(if (antibody == 'fibrinogen' && measure != 'ppositive' && dsname == 'Jamboree' &&  gain < 2500) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'insufficient response'
  },
  silent = TRUE) 
  
  try(if (antibody == 'P-selectin' && measure != 'ppositive' && dsname == 'Jamboree' &&  gain < 9000) {
    DirtyStatus <<- 'Dirty'
    Error <<- 'insufficient response'
  },
  silent = TRUE) 
  
    return(c(DirtyStatus, Error))
}
