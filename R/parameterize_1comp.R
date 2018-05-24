
parameterize_1comp <- function(chem.cas=NULL,chem.name=NULL,species='Human',default.to.human=F, clint.values=NULL,
                               override_httk_param = NULL) {
  
 physiology.data <- physiology.data
   
  if(is.null(chem.cas) & is.null(chem.name)) stop('Must specify chem.name or chem.cas')
  params <- list()
  params[['Vdist']] <- calc_vdist(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human,suppress.messages=T)
  params[['kelim']] <- calc_elimination_rate(chem.cas=chem.cas,chem.name=chem.name,species=species,suppress.messages=T,default.to.human=default.to.human, clint.values=clint.values)
  params[['kgutabs']] <- 1
  params[['Rblood2plasma']] <- calc_rblood2plasma(chem.cas=chem.cas,chem.name=chem.name,species=species,default.to.human=default.to.human)
  params[['million.cells.per.gliver']] <- 110
  
  
  # Check the species argument for capitilization problems and whether or not it is in the table:  
  if (!(species %in% colnames(physiology.data)))
  {
    if (toupper(species) %in% toupper(colnames(physiology.data)))
    {
      phys.species <- colnames(physiology.data)[toupper(colnames(physiology.data))==toupper(species)]
    } else stop(paste("Physiological PK data for",species,"not found."))
  } else phys.species <- species
  
  # Load the physiological parameters for this species
  this.phys.data <- physiology.data[,phys.species]
  names(this.phys.data) <- physiology.data[,1]
  
  # LASER addition: Override physiological parameters
  if(!is.null(override_httk_param)) {
    params_to_update <- names(this.phys.data)[names(this.phys.data) %in% names(override_httk_param)]
    if(length(params_to_update) > 0) {
      for(cpar in params_to_update) {
        if(is.null(override_httk_param[[cpar]]["mean"]) || is.na(override_httk_param[[cpar]]["mean"]))
          this.phys.data[cpar] <- lognormal_var(this.phys.data[cpar], 
                                                override_httk_param[[cpar]]["cv"])
        else
          this.phys.data[cpar] <- lognormal_var(override_httk_param[[cpar]]["mean"], 
                                                override_httk_param[[cpar]]["cv"])
      }
    }
  }
  
  params[['hematocrit']] <- this.phys.data[["Hematocrit"]]
  
  if(is.null(chem.cas)) chem.cas <- get_chem_id(chem.name=chem.name)[['chem.cas']]
  params[['MW']] <- get_physchem_param("MW",chem.CAS=chem.cas)
  
  Fgutabs <- try(get_invitroPK_param("Fgutabs",species,chem.CAS=chem.cas),silent=T)
  if (class(Fgutabs) == "try-error") Fgutabs <- 1
  
  params[['Fgutabs']] <- Fgutabs

  return(params)
}
