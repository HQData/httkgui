solve_mixtures <- function(parametersA, parametersB, parametersint,
                           method="lsoda",rtol=1e-8,atol=1e-12, 
                           use.amounts = F,
                           daily.dose = c(1, 1), dose = NULL, doses.per.day = NULL,
                           ...) {
  if (use.amounts)
    CompartmentsToInitialize <-c("Agutlumen","Aart","Aven","Alung","Agut","Aliver","Akidney","Arest")
  else
    CompartmentsToInitialize <-c("Agutlumen","Cart","Cven","Clung","Cgut","Cliver","Ckidney","Crest")
  
  inputA <- solve_pbtk(parameters = parametersA, prepare_only = TRUE, 
                       daily.dose = daily.dose[1], dose = dose[1], doses.per.day = doses.per.day[1], ...)
  inputB <- solve_pbtk(parameters = parametersB, prepare_only = TRUE, 
                       daily.dose = daily.dose[2], dose = dose[2], doses.per.day = doses.per.day[2], ...)
  names(inputA$y) <- paste0("C1_", names(inputA$y))
  names(inputB$y) <- paste0("C2_", names(inputB$y))
  Outputs <- c("Cgut", "Cliver", "Cven", "Clung", "Cart", "Crest", "Ckidney", "Cserum", "Aserum")
  out <- deSolve::ode(
             y = c(inputA$y, inputB$y),
             times = inputA$times,
             func  = "derivs_mixture",
             parms = c(inputA$parms, inputB$parms, parametersint),
             method=method,rtol=rtol,atol=atol,dllname="httkgui",initfunc="initmod_mixture",
             nout=18, outnames = c(paste0("C1_", Outputs), paste0("C2_", Outputs)))
  

  
  # class(out) <- c('matrix','deSolve')
  
  colnames(out)[[which(colnames(out)=='C1_Aserum')]] <- 'C1_Aplasma'
  colnames(out)[[which(colnames(out)=='C2_Aserum')]] <- 'C2_Aplasma'
  colnames(out)[[which(colnames(out)=='C1_Cserum')]] <- 'C1_Cplasma'
  colnames(out)[[which(colnames(out)=='C2_Cserum')]] <- 'C2_Cplasma'
  
  if (use.amounts)
  {
    CompartmentsToInitialize <-c("Agutlumen","Aart","Aven","Alung","Agut","Aliver","Akidney","Arest", "Aplasma")
  } else {
    CompartmentsToInitialize <-c("Agutlumen","Cart","Cven","Clung","Cgut","Cliver","Ckidney","Crest", "Cplasma")
  }
  
  out1 <- out[,c("time", paste0("C1_", c(CompartmentsToInitialize,"Ametabolized","Atubules", "AUC")))]
  out2 <- out[,c("time", paste0("C2_", c(CompartmentsToInitialize,"Ametabolized","Atubules", "AUC")))]
  colnames(out1) <- sapply(colnames(out1), function(x) {ifelse(x != "time", substr(x, 4, nchar(x)), x)})
  colnames(out2) <- sapply(colnames(out2), function(x) {ifelse(x != "time", substr(x, 4, nchar(x)), x)})
  return(list(out1, out2))
}
