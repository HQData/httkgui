#log-normal variation function according to CV
lognormal_var <- function(x, cv) {
  if(x == 0) return(x)
  if(cv == 0) return(x)
  xsd <- sqrt(log(cv^2 + 1))
  rlnorm(1, log(x) - (xsd^2)/2, xsd)
}

# inits a list of parameters
# param_to_vary is a string vector of names
# multiplier/cv single number, vector or named vector
vary_inits_once <- function(inits, param_to_vary, multiplier = 1, cv = 0) {
  param_to_vary <- as.vector(param_to_vary) #to ensure we're not dealing with factor
  if(!all(param_to_vary %in% names(inits)))
    stop("Unable to change parameter values: not all parameters to vary are present in supplied inits vector")
  value_mean <- unlist(inits[param_to_vary])
  solve_list <- list()
  
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

# generate subpopulation for MC simulations
# - 'base' mean values, some of them can have CV defined
# - for now I assume that multiplier is not defined here but before passing inits
# - CV can be applied before or after parameterization of the model
# individual subpopulations (corresponding to e.g. phenotypes of metabolic pathway) are merged later, using another function
generate_subpopulation <- function(param_to_override = NULL, param_to_vary_before = FALSE, param_to_vary_after = data.frame(), N = 100, name = "") {
  # browser()
  # 0/ set up objects and function arguments:
  pbtk_mclist <- list("inits"=list(), "pbtk_result"=list(), "halflife"=list(), "AUC"=list(), "Cmax"=list())
  param_list <- list("chem.cas"=NULL,"chem.name" = input$compound, 
                     "species" = input$species, "default.to.human" = F,
                     "tissuelist" = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                     "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05, 
                     "monte.carlo" = FALSE, "monte.carlo.cv" = mc_cv(), "monte.carlo.log" = input$mc_use_log,
                     "override_httk_param" = param_to_override)
  
  # variation on CL done this way is now obsolete: we apply it later
  # if(input$cv.clh > 0) 
  #   param_list$clh.cv <- input$cv.clh
  
  if(param_to_vary_before)
    param_list$monte.carlo <- TRUE
  
  if(input$use_cas) {
    param_list$chem.cas <- input$cas
    param_list$chem.name <- NULL
  }
  
  solve_list <- list(
    chem.name=input$compound,
    plots=FALSE, suppress.messages = TRUE,
    output.units = input$solve.output.units, iv.dose = input$solve.iv.dose,
    tsteps = input$solve.tsteps, days = input$solve.days
  )
  if(input$dose_type == "daily dose")
    solve_list$daily.dose=input$solve.daily.dose
  if(input$dose_type == "per dose + doses/day") {
    solve_list$doses.per.day=input$solve.doses
    solve_list$dose=input$solve.dose
  }
  
  # which_are_inits <- custom_param_values$parameter %in% names(parameter_names)
  # which_are_additional <- custom_param_values$parameter %in% names(additional_parameters)
  withProgress(message = "Generating results", min = 0, max = 1, {
    for(i in 1:N) {
      # 1/vary parameters when parameterizing PBTK:
      # --- WIP ---
      
      # 2/parameterize
      inits <- do.call(parameterize_pbtk, param_list)
      
      # 3/ vary parameters after parameterization:
      if(nrow(param_to_vary_after) > 0)
        inits <- vary_inits_once(inits, param_to_vary_after$names, 
                                 cv = param_to_vary_after$cv, 
                                 multiplier = param_to_vary_after$multiplier)
      
      # 4/ solve
      solve_list$parameters <- inits
      pbtk_mclist[["inits"]][[i]] <- inits
      pbtk_mclist[["pbtk_result"]][[i]] <- do.call(solve_pbtk, solve_list)
      # extra calculations (to update)
      times <- pbtk_mclist[["pbtk_result"]][[i]][,"time"]
      x <- pbtk_mclist[["pbtk_result"]][[i]][,"Cplasma"]
      wmax <- which.max(x)
      wmin <- which(x < (max(x)/2))
      pbtk_mclist[["halflife"]][[i]] <- times[min(wmin[wmin > wmax])] - times[wmax]
      pbtk_mclist[["AUC"]][[i]] <- llTrapAUC(times, x)
      pbtk_mclist[["Cmax"]][[i]] <- max(x)
      
      incProgress(1/N)
    }
  })
  # browser()
  pbtk_mclist[["name"]] <- name
  pbtk_mclist
}
