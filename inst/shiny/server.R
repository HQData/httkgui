library(ggplot2)
library(dplyr)
library(tidyr)
# library(plotly)


shiny::shinyServer(function(input, output, session) {
  source("plot_functions.R", local = TRUE)
  source("calculate_functions.R", local = TRUE)
  parameter_names <- c(
      "BW" = "Body Weight, kg.",
      "Clmetabolismc" = "Hepatic Clearance, L/h/kg BW.",
      "Fgutabs" = "Fraction of the oral dose absorbed, i.e. the fraction of the dose that enters the gutlumen.",
      "Funbound.plasma" = "Fraction of plasma that is not bound.",
      "Fhep.assay.correction" = "The fraction of chemical unbound in hepatocyte assay using the method of Kilford et al. (2008)",
      "hematocrit" = "Percent volume of red blood cells in the blood.",
      "kdermabs" = "Rate that chemical is transferred from the skin to the blood, 1/h.",
      "Kgut2pu" = "Ratio of concentration of chemical in gut tissue to unbound concentration in plasma.",
      "kgutabs" = "Rate that chemical enters the gut from gutlumen, 1/h.",
      "kinhabs" = "Rate that the chemical is transferred from the lungs to the blood, 1/h.",
      "Kkidney2pu" = "Ratio of concentration of chemical in kidney tissue to unbound concentration in plasma.",
      "Kliver2pu" = "Ratio of concentration of chemical in liver tissue to unbound concentration in plasma.",
      "Klung2pu" = "Ratio of concentration of chemical in lung tissue to unbound concentration in plasma.",
      "Krbc2pu" = "Ratio of concentration of chemical in red blood cells to unbound concentration in plasma.",
      "Krest2pu" = "Ratio of concentration of chemical in rest of body tissue to unbound concentration in plasma.",
      "million.cells.per.gliver" = "Millions cells per gram of liver tissue.",
      "MW" = "Molecular Weight, g/mol.",
      "Qcardiacc" = "Cardiac Output, L/h/kg BW^3/4.",
      "Qgfrc" = "Glomerular Filtration Rate, L/h/kg BW^3/4, volume of fluid filtered from kidney and excreted.",
      "Qgutf" = "Fraction of cardiac output flowing to the gut.",
      "Qkidneyf" = "Fraction of cardiac output flowing to the kidneys.",
      "Qliverf" = "Fraction of cardiac output flowing to the liver.",
      "Rblood2plasma" = "The ratio of the concentration of the chemical in the blood to the concentration in the plasma.",
      "Vartc" = "Volume of the arteries per kg body weight, L/kg BW.",
      "Vgutc" = "Volume of the gut per kg body weight, L/kg BW.",
      "Vkidneyc" = "Volume of the kidneys per kg body weight, L/kg BW.",
      "Vliverc" = "Volume of the liver per kg body weight, L/kg BW.",
      "Vlungc" = "Volume of the lungs per kg body weight, L/kg BW.",
      "Vrestc" = "Volume of the rest of the body per kg body weight, L/kg BW.",
      "Vvenc" = "Volume of the veins per kg body weight, L/kg BW.",
      "Vmax" = "Maximal velocity, []",
      "km" = "Michaelis constant"
  )
  
  additional_parameters <- c(
    "KTS" = "KTS",
    "FR" = "FR",
    "Clint" = "Clint"
  )
  
  
  
  
  observeEvent(input$use_add, {
      updateTabsetPanel(session, "main_panel",
                        selected = ifelse(input$use_add == 1, "add compound", "inputs summary")
      )
  })
  
  # compound information (define population) --------------------------------
  compound_summary <- reactive({
    filter(chem.physical_and_invitro.data, Compound == input$compound) %>%
      select(Compound, CAS, SMILES.desalt, logP, pKa_Donor, pKa_Accept)
  })
  
  output$compound_table <- renderTable({
    compound_summary()
  })
  
  # observeEvent(input$custom_params, {
  #   if(!input$custom_params && exists("custom_param_values"))
  #     rm(custom_param_values, envir=.GlobalEnv)
  # })
  observeEvent(input$population_new_submit, {
    #custom_subpopulation is used just for user display
    #populations_list is created (for actual calculations) later on
    
    #on first click set them to empty!
    if(input$population_new_submit == 1) {
      custom_subpopulation <<- data.frame()
      populations_list <<- list()
    }
    
    #update the table that the user sees
    custom_subpopulation_newdata <- data.frame(
                          "name"=input$population_new_name, 
                          "N"=input$population_new_N, 
                          "type"=input$population_new_vartype, 
                          "multiplier"=input$population_new_multiplier, 
                          "CV"=input$population_new_cv)
    custom_subpopulation <<- rbind(custom_subpopulation, custom_subpopulation_newdata)
    
    #update the list that guides the simulations
    newlist <- list(
      #fill based on the CV and means provided via custom parameters
      param_to_override = list(
        # "Average BW" = c("mean" = 75, "cv" = 0) 
      ),
      param_to_vary_before = FALSE,
      param_to_vary_after = data.frame(
        "names" = c("Clmetabolismc", "CLmetabolism_gut", "CLmetabolism_kidney"),
        "cv" = input$population_new_cv,
        "multiplier" = input$population_new_multiplier), 
      N = input$population_new_N,
      "name"=input$population_new_name)
    if(input$population_new_vartype == "tk_physbio")
      newlist$param_to_vary_before <- TRUE
    
    populations_list[[length(populations_list) + 1]] <<- newlist
    
    #clean the inputs
    updateTextInput(session, "population_new_name", value = paste("Population", length(populations_list) + 1))
    updateNumericInput(session, "population_new_cv", value = .3)
    updateNumericInput(session, "population_new_N", value = 100)
    # updateNumericInput(session, "population_new_vartype", value = 0)
    updateNumericInput(session, "population_new_multiplier", value = 1)
  })
  # custom_param_values <- data.frame("parameter"=c(), "description"=c(), "value"=c(), 
  # "MC 2.5%"=c(), "MC mean"=c(), "MC 97.5%"=c())
  output$custom_subpopulation_table <- renderTable({
    if((input$population_new_submit > 0) && exists("custom_subpopulation"))
      return(custom_subpopulation)
  })
  
  #experiemental data display:
  output$experimental_data_table <- renderTable({
    experimental_data()
  })
  
  output$model_visual <- renderImage({
    list(src = "pbtk_model.png", alt = "PBTK model schematic", width = 300)
  }, deleteFile = FALSE)
  
  # model parameters  ---------
  
  observeEvent(input$add_submit, {
      if(input$use_add) {
            # browser()
        my.new.data <- data.frame(
            'Compound' = input$add_compound, 
            'CAS' = input$add_cas, 
            'MW' = input$add_mw, 
            'logp'= input$add_logp, 
            'funbound' = input$add_funbound, 
            'fgutabs' = input$add_fgutabs, 
            'clint' = input$add_clint,
            'KTS' = input$add_kts, 
            'FR' = input$add_fr, 
            'vmax' = input$add_vmax, 
            'km' = input$add_km, 
            'pKa_donor' = input$add_pka_donor, 
            'pKa_accept' = input$add_pka_accept)
        
        nna.list <- as.list(na.omit(c(
            'Compound' = 'Compound',
            'CAS' = "CAS",
            'MW' = ifelse(!input$add_mw_na, "MW", NA),
            'logP' = ifelse(!input$add_logp_na, "logp", NA),
            'Funbound.plasma' = ifelse(!input$add_funbound_na, "funbound", NA),
            'Fgutabs' = ifelse(!input$add_fgutabs_na, "fgutabs", NA),
            'Clint' = ifelse(!input$add_clint_na, "clint", NA),
            'KTS' = ifelse(!input$add_kts_na, "KTS", NA),
            'FR' = ifelse(!input$add_fr_na, "FR", NA),
            'Vmax' = ifelse(!input$add_vmax_na, "vmax", NA),
            'km' = ifelse(!input$add_km_na, "km", NA),
            'pKa_Donor' = ifelse(!input$add_pka_donor_na, "pKa_donor", NA),
            'pKa_Accept' = ifelse(!input$add_pka_accept_na, "pKa_accept", NA)
            )))
       
        # 'logMA', 'Clint', 'Clint.pValue', 'Funbound.plasma', 'Fgutabs'
        chem.physical_and_invitro.data_new <<- add_chemtable(my.new.data,
                      current.table=chem.physical_and_invitro.data,
                      data.list=nna.list,
                      species=input$species,
                      reference=input$add_reference, overwrite = TRUE)
      } else {
          
      }
      
  })
  
  observeEvent(input$custom_params, {
    if(!input$custom_params && exists("custom_param_values"))
        # browser()
        rm(custom_param_values, envir=.GlobalEnv)
  })
  observeEvent(input$cparams_submit, {
          #add a row:
          newdata <- data.frame("parameter"=input$cparams_select, 
                                "description"=c(parameter_names, additional_parameters)[input$cparams_select], 
                                "value"=input$cparams_value, 
                                "mc.cv"=input$cparams_cv)
          # browser
          if(exists("custom_param_values")) {
              #remove the last existing value if it was in there
              custom_param_values <<- 
                custom_param_values[custom_param_values$parameter != input$cparams_select,]
              custom_param_values <<- rbind(custom_param_values, newdata)
          } else {
              custom_param_values <<- newdata
          }
      
      #clean the inputs
      updateNumericInput(session, "cparams_value", value = 0)
  })
  # custom_param_values <- data.frame("parameter"=c(), "description"=c(), "value"=c(), 
                                    # "MC 2.5%"=c(), "MC mean"=c(), "MC 97.5%"=c())
  output$custom_param_table <- renderDataTable({
      input$cparams_submit
      if(exists("custom_param_values"))
          return(custom_param_values)
  })
  
  mc_cv <- reactive(c(`Total Body Water` = input$cv.water,
                      `Plasma Volume` = input$cv.plasma,
                      `Cardiac Output` = input$cv.cardiac,
                      `Average BW` = input$cv.bw,
                      `Total Plasma Protein` = input$cv.tpp, 
                      `Plasma albumin` = input$cv.albumin,
                      `Plasma a-1-AGP` = input$cv.a1agp,
                      Hematocrit = input$cv.hematocrit,
                      Urine = input$cv.urine,
                      Bile = input$cv.bile,
                      GFR = input$cv.gfr,
                      `Average Body Temperature` = input$cv.abt
                      ))
  
  #this returns one set of parameters
  parameters <- reactive({
    param_list <- list("chem.cas"=NULL,"chem.name" = input$compound, "species" = input$species, 
                       "default.to.human" = F,
                       "tissuelist" = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                       "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05, monte.carlo=FALSE
                       )
    if(input$use_cas) {
        param_list$chem.cas <- input$cas
        param_list$chem.name <- NULL
    }
    if(input$use_add && input$add_submit) {
        chem.physical_and_invitro.data <<- chem.physical_and_invitro.data_new
        param_list$chem.name <- paste(toupper(substr(input$add_compound, 1, 1)), 
                                      substr(input$add_compound, 2, nchar(input$add_compound)), sep="")
    }
        
    input$cparams_submit
    
    inits <- do.call(parameterize_pbtk, param_list)
    
    #update if user supplied custom values
    if(exists("custom_param_values") && nrow(custom_param_values) > 0) {
      #this part deals with values to update in inits ONLY:  
      which_are_inits <- custom_param_values$parameter %in% names(parameter_names)
      which_are_additional <- custom_param_values$parameter %in% names(additional_parameters)
      if(any(which_are_additional)) {
        #forcing of FR, KTS, Clint
        torep <- custom_param_values$value[which_are_additional]
        names(torep) <- custom_param_values$parameter[which_are_additional]
        param_list$override.input <- torep
        inits <- do.call(parameterize_pbtk, param_list) #overwrite previous calc
      }
        
      torep <- custom_param_values$value[which_are_inits]
      names(torep) <- custom_param_values$parameter[which_are_inits]
      inits[names(torep)] <- torep
    }
    return(inits)
  })
  
  #this will depend on results() object (so generate_subpopulation) as that's where inits live
  parameters_summary <- reactive({
    ww <- parameters()
    if(input$output_type == "single" || (input$run == 0)) {
      parameter_df <- data.frame("parameter"=as.character(names(ww)), 
                                 "description" = parameter_names[names(ww)], 
                                 "value"=unlist(ww))
    }
    else if(input$output_type == "mc") {
      if(is.null(results()))
        return(NULL)
      
      parameter_df <- data.frame("parameter"=names(ww), 
                                 "description" = parameter_names[names(ww)], 
                                 stringsAsFactors = F)
      
      #expand the data frame with uncertainty info:
      df <- lapply(results(), function(x) lapply(x[["inits"]], unlist) %>% do.call(rbind, .)) %>% do.call(rbind, .) %>%
        apply(2, function(x) {c("mean"=mean(x, na.rm=T), 
                                "lci"=quantile(x,.025, na.rm=T), 
                                "uci"=quantile(x,.975, na.rm=T))}) %>%
        t() %>% as.data.frame()
      
      parameter_df[["MC 2.5%"]]  <- df[parameter_df$parameter,2]
      parameter_df[["MC mean"]]  <- df[parameter_df$parameter,1]
      parameter_df[["MC 97.5%"]] <- df[parameter_df$parameter,3]
    }
    return(parameter_df)
  })
  
  output$parameters_df <- renderTable({
    parameters_summary()
  })
  
  # calculation of results(single, monte carlo + summary df for MC) ---------------------------------------
  
  
  #results stored in a single reactive object: both MC and a single simulation
  results <- reactive({
    if(input$output_type == "mc") { #mock eventReactive on input$run
      if(input$run == 0)
        return(NULL)
      
      if(input$population_new_submit < 1) {
        showModal(modalDialog(title = "No population defined", 
                              "Please specify at least one group in the Population Variability section"))
        return(NULL)
      }
      
      lapply(populations_list, function(x) do.call(generate_subpopulation, x))
      
    } else if(input$output_type == "single") {
      generate_subpopulation(N = 1)
    }
  })
  
  endpoints <- reactive({
      # ww <- c("Cplasma", paste0("C", input$compartments), "Crest", "Ametabolized", "Atubules", "Agutlumen")
      ww <- c("Cplasma", paste0("C", c("lung", "kidney", "gut", "liver")), "Crest", "Ametabolized", "Atubules", "Agutlumen")
      names(ww) <- ww
      # names(ww) <- c("Plasma", input$compartments, "rest", "metabolized", "tubules", "gut lumen")
      # names(ww) <- c("Plasma", c("lung", "kidney", "gut", "liver"), "rest", "metabolized", "tubules", "gut lumen")
      return(ww)
  })
  
  experimental_data <- reactive({
    inFile <- input$experimental_data_input
    if(is.null(inFile))
      return(NULL)
    
    tab <- read.csv(inFile$datapath)
    if(is.null(tab$variable))
      tab$variable <- "Cplasma"
    
    tab
  })
  
  results_mc_df_v2 <- reactive({
    if(input$output_type == "single") 
      return(NULL)
    if(is.null(results()))
      return(NULL)
      
    lci_value <- (1-input$display_ci)/2
    uci_value <- 1 - (1-input$display_ci)/2
    res <- results()
    withProgress(message = paste0("Calculating mean parameter values together with ", 100*input$display_ci, "% intervals"), {
      lapply(res, function(current_subpop) {
        tab <- do.call(rbind, current_subpop$pbtk_result)
          tab <- tab %>% as.data.frame() %>% 
            gather(time) %>% setNames(c("time", "variable", "value")) %>% #melt variables
            group_by(time, variable) %>% #group to calculate values for everything
            summarise(mean = mean(value), lci = quantile(value, lci_value), uci = quantile(value, uci_value)) %>%
            mutate(name = current_subpop$name)
      })
    })
  })

  # presentation of results -------------------------------------------------
  
  output$choose_plot_ui <- renderUI({
      selectInput("choose_plot", "Choose parameter to plot", endpoints())
  })
  
  output$choose_plot_type_ui <- renderUI({
    input$run #refresh when we run!
    if(exists("populations_list"))
      if((length(populations_list) > 1) && (input$run > 0))
        return(selectInput("choose_plot_type", "Type of display for subpopulations", c("Color different populations" = "group", 
                                                                                       "Facet (separate panels)" = "facet",
                                                                                       "Both" = "both")))
  })
  
  output$results_plot_single <- renderPlot({
      if(!is.null(input$choose_plot)) {
          if(input$output_type == "mc") {
            if(is.null(results_mc_df_v2()) || is.null(input$choose_plot))
              return(NULL)
            
            #display options:
            fvar <- F; gvar <- F
            if(!is.null(input$choose_plot_type)) {
              if(input$choose_plot_type == "group")
                gvar <- T
              if(input$choose_plot_type == "facet")
                fvar <- T
              if(input$choose_plot_type == "both") {
                fvar <- T; gvar <- T }
            }
            tab <- filter(do.call(rbind, results_mc_df_v2()), variable == input$choose_plot)
            return(auto_gg(tab, facet = fvar, grouping = gvar, varname = input$choose_plot, observed = experimental_data()))
          }
          if(input$output_type == "single") {
            # res <- results_single()[["pbtk_result"]][[1]]
            res <- results()[["pbtk_result"]][[1]]
            cd <- which(colnames(res) == input$choose_plot)
            tab <- res[,c(1, cd)] %>% as.data.frame() %>% setNames(c("time", "mean"))
            return(auto_gg(tab, facet = F, grouping = F, varname = input$choose_plot, observed = experimental_data()))
          }
      }
  })
  
  output$validation_results <- renderUI({
    if(!is.null(experimental_data()))
      return(list(
        h3("Validation against experimental data"),
        plotOutput("results_plot_obspred")
      ))
    
    return(NULL)
    
  })
  
  output$results_plot_obspred <- renderPlot({
    if(!is.null(experimental_data())) {
      if(!is.null(input$choose_plot)) {
        if(input$output_type == "mc") {
          tab <- filter(do.call(rbind, results_mc_df_v2()), variable == input$choose_plot)
          obs <- filter(experimental_data(), variable == input$choose_plot)
          if(nrow(obs) == 0)
            return(NULL)
          return(auto_obspred(prediction = tab, observed = obs))
        }
        if(input$output_type == "single") { 
          # res <- results_single()[["pbtk_result"]][[1]]
          res <- results()[["pbtk_result"]][[1]]
          cd <- which(colnames(res) == input$choose_plot)
          tab <- res[,c(1, cd)] %>% as.data.frame() %>% setNames(c("time", "mean"))
          return(auto_obspred(prediction = tab, observed = experimental_data()))
        }
      }
    }
  })

# left panel observers ----------------------------------------------------

  output$results_plot_ui <- renderUI({
      if(input$output_type == "mc")
          plotOutput("results_plot")
      if(input$output_type == "single")
          plotOutput("results_plot", height=200, width= 600)
  })
  
  
  results_numerical_df <- reactive({
    if(input$output_type == "mc") {
      tab <- do.call(rbind, lapply(results(), function(cpop) {
        lci <- (1-input$display_ci)/2
        uci <- 1 - (1-input$display_ci)/2
        
        df <- apply(data.frame(unlist(cpop[["halflife"]]), 
                               unlist(cpop[["AUC"]]), 
                               unlist(cpop[["Cmax"]])), 
                    2, function(x) {
                      c("lci"=quantile(x,lci, na.rm=T), "mean"=mean(x, na.rm=T), "uci"=quantile(x,uci, na.rm=T)) 
                    })    
        # browser()
        rownames(df) <- c(paste0(100*lci, "%"), "mean", paste0(100*uci, "%"))
        colnames(df) <- c("Half-life (Cplasma)", "AUC (Cplasma)", "Cmax (Cplasma)")
        df <- as.data.frame(t(df))
        if(length(custom_subpopulation) > 1)
          df[["model"]] <- cpop[["name"]]
        
        return(df)
      }))
      
      return(tab)
      
    }
    if(input$output_type == "single") {
      # if(!is.null(results_single())) {
      if(!is.null(results())) {
        # times <- results_single()[,"time"]
        # x <- results_single()[,"Cplasma"]
        # wmax <- which.max(x)
        # wmin <- which(x < (max(x)/2))
        # tt <- times[min(wmin[wmin > wmax])] - times[wmax]
        # return(data.frame("Cplasma half-life" = tt,
        #                   "Cplasma Cmax" = max(x),
        #                   "Cplasma AUC" = llTrapAUC(times, x)))
        # return(data.frame("Cplasma half-life" = results_single()[["halflife"]][[1]],
        #                   "Cplasma Cmax" = results_single()[["Cmax"]][[1]],
        #                   "Cplasma AUC" = results_single()[["AUC"]][[1]]))
        return(data.frame("Cplasma half-life" = results()[["halflife"]][[1]],
                          "Cplasma Cmax" = results()[["Cmax"]][[1]],
                          "Cplasma AUC" = results()[["AUC"]][[1]]))
        
      }
    }
  })
  
  output$results_numerical <- renderTable({
    results_numerical_df()
  }, rownames = TRUE, digits = 3)

  
# reporting -----
  
  output$fileDownload <- downloadHandler(
      filename = function() {
          paste("data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
          if(input$output_type == "single")
              # data <- results_single()
              data <- results()
          if(input$output_type == "mc")
              data <- results_mc_df()["mean",,]
          # browser()
          write.csv(data, 
                    file)
      },
      contentType='text/csv'
  )
  
  output$report <- downloadHandler(
    # filename = "report.html",
    # filename = paste0("report.", input$report_format),
    filename = function() {
      paste('tkplate_report', sep = '.', switch(
        input$report_format, PDF = 'pdf', HTML = 'html', Word = 'docx'
      ))
    },
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      # if(input$report_format == "html"){
      #   tempReport <- file.path(tempdir(), "tkplate_report.Rmd")
      #   file.copy("tkplate_report.Rmd", tempReport, overwrite = TRUE)}
      # if(input$report_format == "pdf"){
      #   tempReport <- file.path(tempdir(), "tkplate_report_pdf.Rmd")
      #   file.copy("tkplate_report_pdf.Rmd", tempReport, overwrite = TRUE)}
      # if(input$report_format == "docx"){
      #   tempReport <- file.path(tempdir(), "tkplate_report_docx.Rmd")
      #   file.copy("tkplate_report_docx.Rmd", tempReport, overwrite = TRUE)}
        
      
      src <- normalizePath('tkplate_report.Rmd')
      # temporarily switch to the temp dir, in case you do not have write
      # permission to the current working directory
      owd <- setwd(tempdir())
      on.exit(setwd(owd))
      file.copy(src, 'tkplate_report.Rmd', overwrite = TRUE)
      
      # Set up parameters to pass to Rmd document
      params <- list(name = input$compound,
                     metabolic_route = input$compound_pathway,
                     compound_chars = compound_summary(),
                     pbtk_parameters = parameters_summary(),
                     population_variability = custom_subpopulation,
                     results = results_numerical_df(),
                     plot = filter(do.call(rbind, results_mc_df_v2()), variable == input$choose_plot) %>% 
                       auto_gg(facet = F, grouping = T, varname = input$choose_plot)
      )
      
      
      if(input$dose_type == "daily dose")
        params$dose <- paste("Single dose (mg/kg BW) of", input$solve.daily.dose)
      if(input$dose_type == "per dose + doses/day")
        params$dose <- paste(input$solve.doses.per.day, "doses per day;", input$solve.dose, "per dose/day (mg/kg BW)")
      
      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      # rmarkdown::render(tempReport, output_file = file,
      #                   params = params,
      #                   envir = new.env(parent = globalenv())
      # )
      out <- rmarkdown::render('tkplate_report.Rmd', output_file = file, params = params, envir = new.env(parent = globalenv()),
        output_format = switch(
        input$report_format,
        PDF = rmarkdown::pdf_document(), HTML = rmarkdown::html_document(), Word = rmarkdown::word_document()
      ))
      file.rename(out, file)
    }
  )

  
  
})


