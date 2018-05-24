#summarise_parameters

summarise_parameters <- function(cpop, variable, conf.int = .95, result = "result") {
  if(any(class(cpop) == "matrix"))
    cpop <- list("result" = list(cpop))
  
  variable #stop if missing!
  
  lci <- (1 - conf.int)/2
  uci <- 1 - (1 - conf.int)/2
  ll <- lapply(cpop[[result]], function(cres) {
    times <- cres[, "time"]
    x     <- cres[, variable]
    hf    <- calculate_halflife(times, x)
    AUC   <- llTrapAUC(times, x)
    Cmax  <- max(x)
    Tmax  <- times[which.max(x)]
    return(data.frame(halflife = hf, AUC = AUC, Cmax = Cmax, 
                      Tmax = Tmax))
  })
  df <- apply(bind_rows(ll), 2, function(x) {
    c(lci = quantile(x, lci, na.rm = T), mean = mean(x, na.rm = T), 
      uci = quantile(x, uci, na.rm = T))
  })
  rownames(df) <- c(paste0(100 * lci, "%"), "mean", paste0(100 * 
                                                             uci, "%"))
  colnames(df) <- paste0(c("Half-life (", "AUC (", "Cmax (", 
                           "Tmax ("), variable, ")")
  df <- as.data.frame(t(df))
  df[["model"]] <- cpop[["name"]]
  return(df)
}