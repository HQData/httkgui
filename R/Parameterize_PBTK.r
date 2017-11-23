# This function parameterizes a PBPK model. The argument tissuelist allows the specific tissues parameerized to be customized.
# All tissues not specified by tissuelist are lumped into a rest of body compartment ("Rest")

parameterize_pbtk <- function(chem.cas = NULL,
                              chem.name = NULL,
                              species = "Human",
                              default.to.human = F,
                              tissuelist=list(liver=c("liver"),kidney=c("kidney"),lung=c("lung"),gut=c("gut")),
                              force.human.clint.fub = F,
                              clint.data = NULL,
                              clint.pvalue.threshold = 0.05,
                              override_httk_param = NULL) {
    
  physiology.data <- physiology.data
  # Look up the chemical name/CAS, depending on what was provide:
  out       <- get_chem_id(chem.cas=chem.cas,chem.name=chem.name)
  chem.cas  <- out$chem.cas
  chem.name <- out$chem.name
   
  if(class(tissuelist)!='list') stop("tissuelist must be a list of vectors.") 
  
  # Clint
  # Clint has units of uL/min/10^6 cells
  
  #LASER addition: _kidney and _gut parameters, Vmax and km method of calculation
  # Clint_kidney <- 0; Clint_gut <- 0; #if they're not provided after this line, 
  #                                    #they will simply be set to 0, and not influence results at all
  
  # if(!is.null(clint.data)) { #this is now the default option for providing Clint_kidney and Clint_gut
  #     Clint <- clint.data[["Vmax"]]*clint.data[["km"]]*clint.data[["ISEF"]]
  #     Clint_kidney <- clint.data[["Vmax_kidney"]]*clint.data[["km_kidney"]]*clint.data[["ISEF_kidney"]]
  #     Clint_gut <- clint.data[["Vmax_gut"]]*clint.data[["km_gut"]]*clint.data[["ISEF_gut"]]
  # } else {
  
  # now try to grab Vmax and km - if they're available, use them to recalculate Clint
  # Vmax <- try(get_invitroPK_param("Vmax", species, chem.CAS=chem.cas), silent=TRUE)
  # km   <- try(get_invitroPK_param("km",   species, chem.CAS=chem.cas), silent=TRUE)
  
  #additional clearance parameters
  Vmax         <- get_override_param("Vmax",         0)
  km           <- get_override_param("km",           0)
  Vmax_kidney  <- get_override_param("Vmax_kidney",  0)
  km_kidney    <- get_override_param("km_kidney",    0)
  Vmax_gut     <- get_override_param("Vmax_gut",     0)
  km_gut       <- get_override_param("km_gut",       0)
  Clint_kidney <- get_override_param("Clint_kidney", 0)
  Clint_gut    <- get_override_param("Clint_gut",    0)
  
  if(km_kidney != 0) { #new calculation (LASER)
    if(Clint_kidney != 0)
      warning("Clint_kidney value provided but not used - calculated as Vmax/km instead")
    Clint_kidney <- Vmax_kidney/km_kidney
  }
  if(km_gut != 0) { #new calculation (LASER)
    if(Clint_gut != 0)
      warning("Clint_gut value provided but not used - calculated as Vmax/km instead")
    Clint_gut <- Vmax_gut/km_gut
  }
  
  #code for liver compartment (original httk)
  # if((class(km) != "try-error") && (class(Vmax) != "try-error")) {
  if(km != 0) { #new calculation (LASER)
    if(get_override_param("Clint", 0) != 0)
      warning("Clint value provided but not used - calculated as Vmax/km instead")
    Clint <- Vmax/km
  } else { #original behaviour
    Clint <- try(get_invitroPK_param("Clint",species,chem.CAS=chem.cas),silent=T)
    if ((class(Clint) == "try-error" & default.to.human) || force.human.clint.fub) {
      Clint <- try(get_invitroPK_param("Clint","Human",chem.CAS=chem.cas),silent=T)
      warning(paste(species,"coerced to Human for metabolic clerance data."))
    }
    if ((class(Clint) == "try-error") && (km == 0 && km_kidney == 0 && km_gut == 0))
      stop("Missing metabolic clearance data for given species. Set default.to.human to true to substitute human value.")
    # Check that the trend in the CLint assay was significant:
    Clint.pValue <- get_invitroPK_param("Clint.pValue",species,chem.CAS=chem.cas)
    if (!is.na(Clint.pValue) & Clint.pValue > clint.pvalue.threshold) 
      Clint <- 0
  }
      
  # unitless fraction of chemical unbound with plasma
  fub <- try(get_invitroPK_param("Funbound.plasma",species,chem.CAS=chem.cas),silent=T)
  if ((class(fub) == "try-error" & default.to.human) || force.human.clint.fub) 
  {
    fub <- try(get_invitroPK_param("Funbound.plasma","Human",chem.CAS=chem.cas),silent=T)
    warning(paste(species,"coerced to Human for protein binding data."))
  }
  if (class(fub) == "try-error") stop("Missing protein binding data for given species. Set default.to.human to true to substitute human value.")
  if (fub == 0)
  {
    fub <- 0.005
    warning("Fraction unbound = 0, changed to 0.005.")
  }
  
  Fgutabs <- try(get_invitroPK_param("Fgutabs",species,chem.CAS=chem.cas),silent=T)
  if (class(Fgutabs) == "try-error") Fgutabs <- 1
    
  
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
  
  temp <- this.phys.data[['Average Body Temperature']] 
  # Load the physico-chemical properties:  
  MW <- get_physchem_param("MW",chem.CAS=chem.cas) #g/mol
  pKa_Donor <- suppressWarnings(get_physchem_param("pKa_Donor",chem.CAS=chem.cas))
  pKa_Accept <- suppressWarnings(get_physchem_param("pKa_Accept",chem.CAS=chem.cas))
  Pow <- 10^get_physchem_param("logP",chem.CAS=chem.cas)
  MA <- suppressWarnings(10^(get_physchem_param("logMA",chem.CAS=chem.cas)))
  
  # Predict the PCs for all tissues in the tissue.data table:
  parm <- list(Funbound.plasma=fub,Pow=Pow,pKa_Donor=pKa_Donor,pKa_Accept=pKa_Accept,MA=MA,Fprotein.plasma = 75/1000/1.025,plasma.pH=7.4,temperature=temp)
  PCs <- predict_partitioning_schmitt(parameters=parm)
  # Get_lumped_tissues returns a list with the lumped PCs, vols, and flows:
  lumped_params <- lump_tissues(PCs,tissuelist=tissuelist,species=species)

  outlist <- list()
  # Begin flows:
  #mL/min/kgBW converted to L/h/kgBW:
  QGFRc <- this.phys.data[["GFR"]]/1000*60
  Qcardiacc <- this.phys.data["Cardiac Output"]/1000*60 
  flows <- unlist(lumped_params[substr(names(lumped_params),1,1) == 'Q'])
  
  # LASER update: consider renal clearance flow instead
  # FR  <- try(get_invitroPK_param("FR",  species, chem.CAS=chem.cas), silent=TRUE)
  # KTS <- try(get_invitroPK_param("KTS", species, chem.CAS=chem.cas), silent=TRUE)
  KTS <- get_override_param("KTS", NA)
  FR  <- get_override_param("FR",  NA)
  if(!is.na(KTS) && !is.na(FR))
    QGFRc <- fub*QGFRc + (flows[["Qkidneyf"]] - fub*QGFRc)*(1-exp(- (fub*QGFRc * KTS/(flows[["Qkidneyf"]]-QGFRc)) ))*(1 - FR)
  # End LASER update of KTS and FR.
    
  outlist <- c(outlist,c(
    Qcardiacc = as.numeric(Qcardiacc),
    flows[!names(flows) %in% c('Qlungf','Qtotal.liverf')],
    Qliverf = flows[['Qtotal.liverf']] - flows[['Qgutf']],
    Qgfrc   = QGFRc
  ))
  # end flows  
  
                                                      
  # Begin volumes
  # units should be L/kgBW  
  Vartc = this.phys.data["Plasma Volume"]/(1-this.phys.data["Hematocrit"])/2/1000 #L/kgBW
  Vvenc = this.phys.data["Plasma Volume"]/(1-this.phys.data["Hematocrit"])/2/1000 #L/kgBW

  outlist <- c(outlist,
    Vartc = as.numeric(Vartc),
    Vvenc = as.numeric(Vvenc),
    lumped_params[substr(names(lumped_params),1,1) == 'V'],
    lumped_params[substr(names(lumped_params),1,1) == 'K'])
  
  
  # Create the list of parameters:
  BW <- this.phys.data["Average BW"]
  hematocrit = this.phys.data["Hematocrit"]
  outlist <- c(outlist,list(
    BW = as.numeric(BW),
    kgutabs  = 1, # 1/h
    kinhabs  = 1, # 1/h
    kdermabs = 1, # 1/h
    Funbound.plasma = as.numeric(fub),        # unitless fraction
    hematocrit      = as.numeric(hematocrit), # unitless ratio
    MW = MW)) #g/mol
  
  # Correct for unbound fraction of chemical in the hepatocyte intrinsic clearance assay (Kilford et al., 2008)
 outlist <- c(outlist,
              list(Fhep.assay.correction = calc_fu_hep(Pow,pKa_Donor=pKa_Donor,pKa_Accept=pKa_Accept)))

 #parameters needed to calculate CL metabolism (first - case of liver)
  CLh_parameters <- list(
    Clint=Clint, #uL/min/10^6 cells
    Funbound.plasma=fub, # unitless fraction
    Fhep.assay.correction=outlist$Fhep.assay.correction, 
    million.cells.per.gliver= 110, # 10^6 cells/g-liver
    liver.density= 1.05, # g/mL
    Dn=0.17,BW=BW,
    Vliverc=lumped_params$Vliverc, #L/kg
    Qtotal.liverc=(lumped_params$Qtotal.liverc)/1000*60)

  #differences between liver vs kidney & gut are on 4 parameters:
  CLh_parameters_kidney <- CLh_parameters
  CLh_parameters_kidney[["Clint"]] <- Clint_kidney
  CLh_parameters_kidney[["million.cells.per.gliver"]] <- 13.6 #this is in reality MPPGK but we keep parameter label that relates to liver
  CLh_parameters_kidney[["Vliverc"]] <- lumped_params$Vkidneyc
  CLh_parameters_kidney[["liver.density"]] <- 1.05
  
  CLh_parameters_gut <- CLh_parameters
  CLh_parameters_gut[["Clint"]] <- Clint_gut
  CLh_parameters_gut[["million.cells.per.gliver"]] <- 21.1 #MPPGG
  CLh_parameters_gut[["Vliverc"]] <- lumped_params$Vgutc
  CLh_parameters_gut[["liver.density"]] <- 1.05
  
  #calculate using hepatic clearance  function for all 3 compartments:
  CLh_value      <- as.numeric(calc_hepatic_clearance(hepatic.model="unscaled", 
                                                      parameters=CLh_parameters, suppress.messages=T))
  CLkidney_value <- as.numeric(calc_hepatic_clearance(hepatic.model="unscaled", 
                                                      parameters=CLh_parameters_kidney, suppress.messages=T))
  CLgut_value    <- as.numeric(calc_hepatic_clearance(hepatic.model="unscaled", 
                                                      parameters=CLh_parameters_gut, suppress.messages=T))
  
  outlist <- c(outlist,
     list(Clmetabolismc            = CLh_value,
          CLmetabolism_gut         = CLgut_value,
          CLmetabolism_kidney      = CLkidney_value,
          million.cells.per.gliver = 110,
          Fgutabs                  = Fgutabs         #L/h/kg BW
          )
     ) 
  outlist <- c(outlist,Rblood2plasma=as.numeric(1 - hematocrit + hematocrit * PCs[["Krbc2pu"]] * fub))
 
 
  #LASER addition - extra Vmax and km parameters required by the extended C model:
  outlist <- c(outlist, list("Vmax" = Vmax, "km" = km, "Vmax_kidney" = Vmax_kidney, 
                             "km_kidney" = km_kidney, "Vmax_gut" = Vmax_gut, "km_gut" = km_gut))
  return(outlist)
}