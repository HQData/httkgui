#summarise_population

summarise_population <- function(current_pop, lci_value = .025, uci_value = .975, result = "result"){
  if(all(class(current_pop) == "list")){
      tab <- do.call(rbind, current_pop[[result]])
  }
  else if(any(class(current_pop) == "matrix"))
    tab <- current_pop
  tab <- tab %>% as.data.frame() %>% 
    tidyr::gather(variable, value, -time) %>% 
    # setNames(c("time", "variable", "value")) %>% 
    dplyr::group_by(time, variable) %>% 
    dplyr::summarise(mean = mean(value), lci = quantile(value, lci_value), uci = quantile(value, uci_value))
  
  if(all(class(current_pop) == "list"))
    tab$name <- current_pop[["name"]]
  else
    tab$name <- ""
  
  return(tab)
}
