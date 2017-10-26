# generate subpopulation for MC simulations (for easier usage I call it generate_population)
# - 'base' mean values, some of them can have CV defined
# - for now I assume that multiplier is not defined here but before passing inits
# - CV can be applied before or after parameterization of the model
# individual subpopulations (corresponding to e.g. phenotypes of metabolic pathway) are merged later, using another function


generate_population <- function(compound, species, cas = NULL, use.cas = F,
                                model = "pbtk",
                                solve.output.units = 'uM', solve.iv.dose = F, solve.tsteps = 4, solve.days = 1,
                                dose_type = "daily", solve.daily.dose = 1, solve.doses = NULL, solve.dose = NULL,
                                param_to_override = NULL, param_to_vary_before = FALSE, 
                                param_to_vary_after = data.frame(), N = 100, name = "") {
  # 0/ set up objects and function arguments:
  mclist <- list("inits"=list(), "result"=list(), "halflife"=list(), "AUC"=list(), "Cmax"=list())
  # list of arguments to parameterize
  parameterize_inputs <- list("chem.cas"=NULL,"chem.name" = compound, 
                     "species" = species, "default.to.human" = F)
  
  if(model == "pbtk"){
    parameterize_inputs <- append(parameterize_inputs, list(
                     "tissuelist" = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                     "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05, 
                     "monte.carlo" = FALSE, "override_httk_param" = param_to_override))
    #this is legacy functionality, no need to go there for now!
    if(param_to_vary_before)
      parameterize_inputs$monte.carlo <- TRUE
  } else if(model == "1comp") {
    # nothing to add for now, 
    #in future should test custom Clint values
  } else {
    stop("Argument 'model' has to be either 'pbtk' or '1comp'.")
  }
  if(use.cas) {
    parameterize_inputs$chem.cas <- cas
    parameterize_inputs$chem.name <- NULL
  }
  
  # 1/ inputs for the solver:
  # (same for PBTK and 1comp cases)
  solve_list <- list(
    chem.name=compound,
    plots=FALSE, suppress.messages = TRUE,
    output.units = solve.output.units, iv.dose = solve.iv.dose,
    tsteps = solve.tsteps, days = solve.days
  )
  if(dose_type == "daily")
    solve_list$daily.dose=solve.daily.dose
  if(dose_type == "per day") {
    solve_list$doses.per.day=solve.doses
    solve_list$dose=solve.dose
  }
  
  for(i in 1:N) {
    # a/ parameterize
    if(model == "pbtk")
      inits <- do.call(parameterize_pbtk,  parameterize_inputs)
    if(model == "1comp")
      inits <- do.call(parameterize_1comp, parameterize_inputs)
      
    # b/ vary parameters after parameterization:
    if(nrow(param_to_vary_after) > 0)
      inits <- vary_inits(inits, param_to_vary_after$names, 
                          cv = param_to_vary_after$cv, 
                          multiplier = param_to_vary_after$multiplier)
    
    # c/ solve
    solve_list$parameters  <- inits
    mclist[["inits"]][[i]] <- inits
    if(model == "pbtk")
      mclist[["result"]][[i]] <- do.call(solve_pbtk,  solve_list)
    if(model == "1comp")
      mclist[["result"]][[i]] <- do.call(solve_1comp, solve_list)
    
    # d/ extra calculations (to update for AUC, halflife, Cmax of all variables?)
    if(model == "pbtk")
      varname <- "Cplasma"
    if(model == "1comp")
      varname <- "Ccompartment"
    times <- mclist[["result"]][[i]][,"time"]
    x <- mclist[["result"]][[i]][,varname]
    wmax <- which.max(x)
    wmin <- which(x < (max(x)/2))
    mclist[["halflife"]][[i]] <- times[min(wmin[wmin > wmax])] - times[wmax]
    mclist[["AUC"]][[i]] <- llTrapAUC(times, x)
    mclist[["Cmax"]][[i]] <- max(x)
    
    # if(shiny)
      # shiny::incProgress(1/N)
  }
  mclist[["name"]] <- name
  mclist
}
