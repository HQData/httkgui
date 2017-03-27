llTrapAUC <- function(time, values) {
    dt <- time[2:length(time)] - time[1:(length(time) - 1)]
    c2 <- values[2:length(values)]
    c1 <- values[1:(length(values) - 1)]
    if(length(c1) != length(dt))
        stop("In AUC calculation length of times and values do not match.")
    
    # cat(dt)
    ww <- apply(data.frame(dt, c1, c2), 1, function(x) {
      if(x[3] >= x[2])
          return(x[1]*(x[3]+x[2])/2)
      if(x[3] < x[2])
          return(x[1]*(x[2] - x[3])/(log(x[2]) - log(x[3])))
    })
    
    return(sum(ww))
}
