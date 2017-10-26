#log-normal variation function according to CV
lognormal_var <- function(x, cv) {
  if(x == 0) return(x)
  if(cv == 0) return(x)
  xsd <- sqrt(log(cv^2 + 1))
  rlnorm(1, log(x) - (xsd^2)/2, xsd)
}