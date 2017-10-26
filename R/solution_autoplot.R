# visualisation functions -------------------------------------------------
solution_autoplot <- function(ggdata, type="PBTK", nameval = "all", facet = FALSE, grouping = FALSE, varname = "Cplasma", observed = NULL) {
  # if(!is.null(facet)) ggdata$facet <- ggdata[[facet]]
  # named <- paste0("output/", compound_name, "_", type, "_plot_", nameval, ".jpg")
  if(is.null(ggdata$name)) ggdata$name <- "All"
  if(is.null(ggdata$lci) || is.null(ggdata$uci)) {
    ggdata$lci <- ggdata$mean
    ggdata$uci <- ggdata$mean
  }
  gg <- ggplot2::ggplot(ggdata, aes(y = mean, x = time)) + 
  {if(grouping) geom_line(size = 1.5, aes(color=name, group = name)) } +
  {if(!grouping) geom_line(size = 1.5) } +
  {if(grouping) geom_ribbon(aes(min = lci, max=uci, fill=name, group = name), alpha=.1) } +
  {if(grouping) geom_line(aes(y = lci, color=name, group = name), linetype = "dashed") } +
  {if(grouping) geom_line(aes(y = uci, color=name, group = name), linetype = "dashed") } +
  {if(!grouping) geom_ribbon(aes(min = lci, max=uci), alpha=.25) } +
    ylab(varname) + 
    # ggtitle(paste0(varname, " in ", compound_name, " - ", type, " model")) + 
    {if(max(ggdata$time) < 1.5) scale_x_continuous(breaks = seq(0, 1, length = 25), labels = 0:24) } +
    {if(max(ggdata$time) < 1.5) xlab("Time (hours)") } + 
    {if(max(ggdata$time) >= 1.5) xlab("Time (days)") } + 
    theme_bw(base_size = 20) + theme(legend.position = "top") +
    scale_fill_discrete(guide = FALSE) +
    { if(!grouping)  scale_color_discrete(guide = FALSE) } +
    { if(facet) facet_grid(. ~ name) }
  
  #observed is a df with cols: name, time, mean, lower (optional), upper (optional), observed_group (optional)
  if(is.null(observed$lower) || is.null(observed$upper)) {
    observed$lower <- observed$mean
    observed$upper <- observed$mean
  }
  if(!is.null(observed)) {
    if(!is.null(observed$name)) {
      if(!grouping) compare <- observed %>% dplyr::filter(name == "All")
      if(grouping) compare <- observed %>% dplyr::filter(name != "All")
    } else {
      compare <- observed
    }
    if(nrow(compare) > 0) {
      if(is.null(observed$observed_group)) {
        gg <- gg +
          geom_errorbar(data = compare, aes(ymin = lower, ymax = upper, x = time), width = 0) +
          geom_point(data = compare, aes(x = time, y = mean), size = 3)
      } else {
        gg <- gg +
        {if(grouping) geom_point(data = compare, aes(x = time, y = mean, color = name), size = 3) } +
        {if(!grouping) geom_point(data = compare, aes(x = time, y = mean, color = observed_group), size = 3) } +
        {if(grouping) geom_errorbar(data = compare, aes(ymin = lower, ymax = upper, x = time, color = name), width = 0) } +
        {if(!grouping) geom_errorbar(data = compare, aes(ymin = lower, ymax = upper, x = time, color = observed_group), width = 0)}
      }
    }
  }
  
  # ggsave(plot = gg, filename = "test.jpg", units="cm", width=17, height=12)
  # ggsave(plot = gg, filename = named, units="cm", width=17, height=12)
  return(gg)
}
