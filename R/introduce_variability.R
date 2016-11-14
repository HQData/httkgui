
# LASER-added function: introducing variability to the whole matrix of interest
# simple scheme where we use .mean and .cv columns to modify ALL values
# computation-intensive for longer sims
# (of course, this function itself is a few secs for 10,000 runs, but what about solve?)

introduce_variability <- function(input.mean, input.sd=NULL, cv=NULL, dist=NULL, positive=TRUE, input.name.col = "Parameter") {
    #vary all parameters where CV column matching names of means can be found
    whichnum <- lapply(input.mean, class) == "numeric"
    coln <- colnames(input.mean)
    
    if(is.null(input.sd) & is.null(cv))
        stop("Please specify either SD table or CV for all parameters in the input table")
    
    if(!is.null(input.sd)) {
        ii <- 1
        if(!is.null(cv)) message("Both SD and CV specified for introducing variability. Using SD as default.")
        return(
            data.frame(lapply(input.mean, function(x) {
                res <- x
                if(is.numeric(x))
                    res <- rnorm(length(x), x, input.sd[[coln[ii]]])
                ii <<- ii + 1
                return(res)
            }))
        )
    }
    
    if(!is.null(cv)) {
        if(!is.null(names(cv))) {
            if(!all(names(cv) %in% input.mean[[input.name.col]]))
                stop("CV not specified for all parameters in the table.")
            cv <- cv[input.mean[[input.name.col]]]
        }
        return(
            data.frame(lapply(input.mean, function(x) {
                res <- x
                if(is.numeric(x))
                    res <- rnorm(length(x), x, cv*x)
                    res[res < 0] <- 0 #this simple fix mightstop crashes but generally we should consider log dist.
                return(res)
            }))
        )        
    }
}