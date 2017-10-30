calculate_halflife <- function(times, y) {
  wmax <- which.max(y)
  wmin <- which(y < (max(y)/2))
  
  #if nothing meets the halflife condition...
  if(length(wmin[wmin > wmax]) == 0)
    return(NA)
  
  times[min(wmin[wmin > wmax])] - times[wmax]
}