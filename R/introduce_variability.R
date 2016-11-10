
# LASER-added function: introducing variability to the whole matrix of interest
# simple scheme where we use .mean and .cv columns to modify ALL values
# computation-intensive for longer sims
# (of course, this function itself is a few secs for 10,000 runs, but what about solve?)

introduce_variability <- function(input.mean, input.sd, dist=NULL, positive=TRUE) {
    #vary all parameters where CV column matching names of means can be found
    whichnum <- lapply(input.mean, class) == "numeric"
    coln <- colnames(input.mean)
    ii <- 1
    
    data.frame(lapply(input.mean, function(x) {
        res <- x
        if(is.numeric(x))
            res <- rnorm(length(x), x, input.sd[[coln[ii]]])
        ii <<- ii + 1
        return(res)
    }))
}