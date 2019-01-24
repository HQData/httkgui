sensitivity_test <- function(
    param_values,
    param_to_vary,
    method="uniform10",
    #use.amounts from solve_pbtk give us A...
    state = c("Agutlumen"=332.9498,"Aart"=0,"Aven"=0,"Alung"=0,"Agut"=0,"Aliver"=0,"Akidney"=0,"Arest"=0),
    days=10,
    tsteps = 4,
    times = round(seq(0, days, 1/(24*tsteps)),8),
    Nsim = 1000
) {
    if(is.null(param_to_vary) || (length(param_to_vary) == 0))
      stop("No parameters to vary.")
    
    sensitivity_result <- list()
    for(i in 1:Nsim) {
      if(method=="uniform10")
        param_values[param_to_vary] <- sapply(param_values[param_to_vary], function(x) runif(1, .9*x, 1.1*x))
      # browser()
      
      # state <-initState(param_values,state) # not needed?
      # out <- 1
      out <- ode(y = state,
                 times = times,
                 func="derivs",
                 parms=param_values,
                 method="lsoda",rtol=1e-8,atol=1e-12,
                 dllname="httkgui",initfunc="initmod",
                 nout=length(Outputs),outnames=Outputs)
      class(out) <- c('matrix','deSolve')
      sensitivity_result[[i]] <- out
    }
    
    return(sensitivity_result)
    
    
}