# inits a list of parameters
# param_to_vary is a string vector of names
# multiplier/cv single number, vector or named vector
vary_inits <- function(inits, param_to_vary, mean = NULL, multiplier = 1, cv = 0) {
  param_to_vary <- as.vector(param_to_vary) #to ensure we're not dealing with factor
  if(!all(param_to_vary %in% names(inits)))
    stop("Unable to change parameter values: not all parameters to vary are present in supplied inits vector")
  value_mean <- unlist(inits[param_to_vary])
  solve_list <- list()
  
  if(!is.null(mean)){
    if(is.null(names(mean)))
      names(mean) <- param_to_vary
    #if some NA's have been provided, just override them with defaults
    mean[is.na(mean)] <- value_mean[names(mean)[is.na(mean)]]
    value_mean[names(mean)] <- mean
  }
  
  if(is.null(multiplier))
    multiplier <- 1
  
  #if there's just 1 value, we assume that it's supposed to be the same for all
  if(length(cv) == 1) cv <- rep(cv, length(param_to_vary))
  if(length(multiplier) == 1) multiplier <- rep(multiplier, length(param_to_vary))
  #if no names were supplied for this vector, we assume that order is same as param_to_vary vector
  if(is.null(names(multiplier))) names(multiplier) <- param_to_vary
  if(is.null(names(cv))) names(cv) <- param_to_vary
  
  #apply CV and multiplier
  for(cpar in param_to_vary)
    inits[[cpar]] <- as.numeric(lognormal_var(multiplier[cpar]*value_mean[cpar], cv[cpar]))
  
  return(inits)
}
