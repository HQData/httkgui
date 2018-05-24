# generate (sub)population for MC simulations (for easier usage I call it generate_population)
# - 'base' mean values, some of them can have CV defined
# - for now I assume that multiplier is not defined here but before passing inits
# - CV can be applied before or after parameterization of the model
# individual subpopulations (corresponding to e.g. phenotypes of metabolic pathway) are merged later, using another function


generate_population <- function(compound, species, cas = NULL, use.cas = F,
                                model = "pbtk", default.to.human = T,
                                #all the inputs passed on to model solver:
                                output.units = 'uM', iv.dose = F, tsteps = 4, days = 1,
                                daily.dose = 1, dose = NULL, doses.per.day = NULL,
                                #choice of custom values and other options:
                                param_to_override = NULL, param_to_vary_after = data.frame(), 
                                parametersint = NULL,
                                N = 100, name = "",
                                custom_inits = NULL) {
  # 0/ set up objects and function arguments:
  if(model == "mixture")
    mclist <- list("initsA"=list(), "initsB" = list(), "resultA"=list(), "resultB"=list(), 
                   "parametersint" = list(),
                   "halflife"=list(), "AUC"=list(), "Cmax"=list())
  else 
    mclist <- list("inits"=list(), "result"=list(), "halflife"=list(), "AUC"=list(), "Cmax"=list())
  
  # list of arguments to parameterize
  parameterize_inputs <- list("chem.cas"=NULL,"chem.name" = compound, 
                     "species" = species, "default.to.human" = default.to.human)
  
  if(model == "pbtk"){
    parameterize_inputs <- append(parameterize_inputs, list(
                     "tissuelist" = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                     "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05, "override_httk_param" = param_to_override))
                     
  } else if(model == "1comp") {
    parameterize_inputs <- append(parameterize_inputs, list(
      "override_httk_param" = param_to_override))
  } else if(model == "mixture") {
    parameterize_inputsA <- list("chem.cas"=NULL,"chem.name" = compound[1], 
                                 "species" = species, "default.to.human" = default.to.human,
                                 "tissuelist" = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                                 "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05,
                                 "override_httk_param" = param_to_override[[1]])
    parameterize_inputsB <- list("chem.cas"=NULL,"chem.name" = compound[2], 
                                 "species" = species, "default.to.human" = default.to.human,
                                 "tissuelist" = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                                 "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05,
                                 "override_httk_param" = param_to_override[[2]])
  } else {
    stop("Argument 'model' has to be 'pbtk' or '1comp' or 'mixture'.")
  }
  if(use.cas) {
    parameterize_inputs$chem.cas <- cas
    parameterize_inputs$chem.name <- NULL
  }
  
  # 1/ inputs for the solver:
  # (same for PBTK and 1comp cases)
  solve_list <- list(
    "chem.name"=compound,
    "plots"=FALSE, "suppress.messages" = TRUE,
    "output.units" = output.units, "iv.dose" = iv.dose,
    "tsteps" = tsteps, "days" = days
  )
  #for dose for now no complications: same behaviour as in solve_pbtk/solve_1comp
  solve_list$daily.dose    <- daily.dose
  solve_list$dose          <- dose
  solve_list$doses.per.day <- doses.per.day
  
  for(i in 1:N) {
    # a/ parameterize
    if(model == "pbtk")
      inits <- do.call(parameterize_pbtk,  parameterize_inputs)
    if(model == "1comp")
      inits <- do.call(parameterize_1comp, parameterize_inputs)
    if(model == "mixture"){
      initsA <- do.call(parameterize_pbtk, parameterize_inputsA)
      initsB <- do.call(parameterize_pbtk, parameterize_inputsB)
    }
    
    # b/ vary parameters after parameterization:
    if(model != "mixture") {
      if(!is.null(custom_inits)) {
        solve_list$parameters  <- custom_inits[[i]]
        mclist[["inits"]][[i]] <- custom_inits[[i]]
      } else {
        if(nrow(param_to_vary_after) > 0)
          inits <- vary_inits(inits, 
                              param_to_vary_after$names, 
                              mean       = param_to_vary_after$mean,
                              cv         = param_to_vary_after$cv, 
                              multiplier = param_to_vary_after$multiplier)
        solve_list$parameters  <- inits
        mclist[["inits"]][[i]] <- inits
      }
    } else {
      if(!is.null(custom_inits)) {
        solve_list$parametersA   <- custom_inits[[1]][[i]]
        solve_list$parametersB   <- custom_inits[[2]][[i]]
        mclist[["initsA"]][[i]]  <- custom_inits[[1]][[i]]
        mclist[["initsB"]][[i]]  <- custom_inits[[2]][[i]]
        solve_list$parametersint <- custom_inits[[3]][[i]]
        mclist[["parametersint"]][[i]] <- custom_inits[[3]][[i]]
      } else {
        initsA <- vary_inits(initsA, 
                             param_to_vary_after[[1]]$names, 
                             mean       = param_to_vary_after[[1]]$mean,
                             cv         = param_to_vary_after[[1]]$cv, 
                             multiplier = param_to_vary_after[[1]]$multiplier)
        initsB <- vary_inits(initsB, 
                             param_to_vary_after[[2]]$names, 
                             mean       = param_to_vary_after[[2]]$mean,
                             cv         = param_to_vary_after[[2]]$cv, 
                             multiplier = param_to_vary_after[[2]]$multiplier)
        solve_list$parametersA   <- initsA
        solve_list$parametersB   <- initsB
        mclist[["initsA"]][[i]] <- initsA
        mclist[["initsB"]][[i]] <- initsB
        parametersint <- vary_inits(parametersint, 
                                    param_to_vary_after[["interaction"]]$names, 
                                    mean       = param_to_vary_after[["interaction"]]$mean,
                                    cv         = param_to_vary_after[["interaction"]]$cv, 
                                    multiplier = param_to_vary_after[["interaction"]]$multiplier)
        solve_list$parametersint <- parametersint
        mclist[["parametersint"]][[i]] <- parametersint
      }
      
      
    }
    
    
    # c/ solve
    if(model == "pbtk")
      mclist[["result"]][[i]] <- do.call(solve_pbtk,  solve_list)
    if(model == "1comp")
      mclist[["result"]][[i]] <- do.call(solve_1comp, solve_list)
    if(model == "mixture"){
      res <- do.call(solve_mixtures, solve_list)
      mclist[["resultA"]][[i]] <- res[[1]]
      mclist[["resultB"]][[i]] <- res[[2]]
    }
    # d/ extra calculations (to update for AUC, halflife, Cmax of all variables?)
    # now nothing is precalculated, we do it in summarise_statistics
    
  }
  mclist[["name"]] <- name
  mclist
}


