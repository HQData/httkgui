
#summarise_population
summarise_population <- function(current_pop, lci_value = .025, uci_value = .975){
  tab <- do.call(rbind, current_pop$result)
  tab <- tab %>% as.data.frame() %>% 
    tidyr::gather(time) %>% 
    setNames(c("time", "variable", "value")) %>% #melt variables
    dplyr::group_by(time, variable) %>% #group to calculate values for everything
    dplyr::summarise(mean = mean(value), lci = quantile(value, lci_value), uci = quantile(value, uci_value)) %>%
    dplyr::mutate(name = current_pop$name)
}