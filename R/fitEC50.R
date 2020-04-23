#' Fit a Hill curve to the data
#'
#' @param L, data (x-axis)
#'
#' @param Y, data (y-axis)
#'
#' @param measure, ppositive or MedianFI
#'
#' @return The ec50, the hill coefficient and goodness of fit
fitEC50 <- function(L, Y, measure) {
  bottom <- min(Y)
  top <- max(Y)
  
  suppressWarnings(fitted <-
                     try(fit <- nlsLM(
                       Y ~ bottom + (top - bottom) / (1 + (ec50 / L) ^ hill),
                       start = c(hill = 6, ec50 = 14),
                       lower = c(hill = 0, ec50 = 0),
                       upper = c(hill = 150, ec50 = 100),
                       control = nls.lm.control(maxiter = 100)
                     ),
                     silent = TRUE)
  )
  
  if (class(fitted) != "try-error")
  {
    EC50 <<- coef(fitted)["ec50"]
    Hill <<- coef(fitted)["hill"]
    
    if (measure != 'ppositive')
    {
      Goodness <- sum(residuals(fit) ^ 2)  /  top
    }  else {
      Goodness <- sum(residuals(fit) ^ 2)
    }
  } else {
    EC50 <<- 0
    Hill <<- 0
    Goodness <<- 0
  }
  
  return(c(EC50, Hill, Goodness))
  
}
