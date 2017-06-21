# This function calculats an elimination rate for a one compartment model where 
# eelimination is entirely due to metablism by the liver and glomerular filtration
# in the kidneys.
calc_total_clearance<- function(chem.cas=NULL,chem.name=NULL,parameters=NULL,species="Human",suppress.messages=F,
                                 default.to.human=F, clint.values=NULL)
{
    if(is.null(parameters)) parameters <- parameterize_steadystate(chem.cas=chem.cas, chem.name=chem.name, species=species,default.to.human=default.to.human)
    
    
    Qgfrc <- get_param("Qgfrc",parameters,"calc_Css") / parameters[['BW']]^0.25 #L/h/kgBW
    fub <- get_param("Funbound.plasma",parameters,"calc_Css") # unitless fraction
    
    # CL is now comprised of three components:
    CLparameters <- parameters
    CLh_parameters_kidney <- CLparameters
    CLh_parameters_gut <- CLparameters
    
    if(is.null(clint.values)) {
      Clint_kidney <- 0
      Clint_gut <- 0
    } else {
      Clint_kidney <- clint.values[["Clint_kidney"]]
      Clint_gut <- clint.values[["Clint_gut"]]
      if(!is.null(clint.values[["Clint_liver"]]))
        CLparameters[["Clint"]] <- clint.values[["Clint_liver"]]
    }
      
    CLh_parameters_kidney[["Clint"]] <- Clint_kidney
    CLh_parameters_kidney[["million.cells.per.gliver"]] <- 13.6 #this is in reality MPPGK but we keep parameter label that relates to liver
    CLh_parameters_kidney[["Vliverc"]] <- 0.004190476
    CLh_parameters_kidney[["liver.density"]] <- 1.05
    
    CLh_parameters_gut[["Clint"]] <- Clint_gut
    CLh_parameters_gut[["million.cells.per.gliver"]] <- 3.79 #MPPGG
    CLh_parameters_gut[["Vliverc"]] <- 0.015802847
    CLh_parameters_gut[["liver.density"]] <- 1.05
    
    CLsum <- calc_hepatic_clearance(parameters=CLparameters, suppress.messages=T)
    CLsum <- CLsum + calc_hepatic_clearance(parameters=CLh_parameters_kidney, suppress.messages=T)
    CLsum <- CLsum + calc_hepatic_clearance(parameters=CLh_parameters_gut, suppress.messages=T)
    
    clearance <- Qgfrc*fub + CLsum #L/h/kgBW
    if(!suppress.messages)cat(paste(toupper(substr(species,1,1)),substr(species,2,nchar(species)),sep=''),"total clearance returned in units of L/h/kg BW.\n")
    return(as.numeric(clearance))
}