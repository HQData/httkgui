#for introducing variability in additional parameters that have been added by LASER modifications
get_override_param <- function(param, default = 0) {  
  if(!exists("override_httk_param", parent.frame()))
    return(default)
  else
    override_httk_param <- get("override_httk_param", parent.frame())
  
  if(param %in% names(override_httk_param)){
    cv <- ifelse(!is.null(override_httk_param[[param]]["cv"]), override_httk_param[[param]]["cv"], 0)
    return(as.numeric(lognormal_var(override_httk_param[[param]]["mean"], cv)))
  } else {
    return(default)
  }
}