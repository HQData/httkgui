#log-normal variation function according to CV
lognormal_var <- function(x, cv) {
  if(x == 0) return(x)
  if(cv == 0) return(x)
  xsd <- sqrt(log(cv^2 + 1))
  rlnorm(1, log(x) - (xsd^2)/2, xsd)
}




# this should be incorporated into param_to_vary_before, param_to_vary_after:
# #update if user supplied custom values
# if(exists("custom_param_values") && nrow(custom_param_values) > 0) {
#   #this part deals with values to update in inits ONLY:  
#   which_are_inits <- custom_param_values$parameter %in% names(parameter_names)
#   which_are_additional <- custom_param_values$parameter %in% names(additional_parameters)
#   if(any(which_are_additional)) {
#     #forcing of FR, KTS, Clint
#     torep <- custom_param_values$value[which_are_additional]
#     names(torep) <- custom_param_values$parameter[which_are_additional]
#     param_list$override.input <- torep
#     inits <- do.call(parameterize_pbtk, param_list) #overwrite previous calc
#   }
#   
#   torep <- custom_param_values$value[which_are_inits]
#   names(torep) <- custom_param_values$parameter[which_are_inits]
#   inits[names(torep)] <- torep
# }

