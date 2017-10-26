

#two df's observed and predicted only need mean and time columns... for now
auto_obspred <- function(observed, prediction) {
  df_op <- data.frame()
  for(i in 1:nrow(observed)) {
    t1 <- observed$time[i]
    m1 <- observed$mean[i]
    wm <- which.min(abs(prediction$time - t1))
    # t2 <- predicted$time[wm]
    m2 <- prediction$mean[wm]
    df_op <- rbind(df_op, data.frame("time" = t1, "observed" = m1, "predicted" = m2))
  }
  coordmax <- max(c(observed$mean, prediction$mean))
  gg <- ggplot(df_op) + geom_point(aes(x=observed, y=predicted), shape = 18, size = 2) + 
    geom_abline(intercept = 0, slope = 1) + 
    geom_abline(intercept=0, slope = .1, linetype = 3) +
    geom_abline(intercept=0, slope = .33, linetype = 2) +
    geom_abline(intercept=0, slope = 10, linetype = 3) +
    geom_abline(intercept=0, slope = 3, linetype = 2) +
    coord_cartesian(xlim = c(0, coordmax * 1.1), ylim = c(0, coordmax * 1.1)) +
    theme_minimal(base_size = 20) +
    # xlab(paste("observed", output_var_pbtk, default_output_unit)) +
    # ylab(paste("predicted", output_var_pbtk, default_output_unit)) +
    ggtitle("Observed vs predicted plot for the model") +
    labs(subtitle = "(dotted line = 10-fold difference, dashed = 3-fold)")
  
  return(gg)
}
