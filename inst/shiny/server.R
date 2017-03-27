

shiny::shinyServer(function(input, output, session) {

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
    
    observeEvent(input$use_add, {
        updateTabsetPanel(session, "main_panel",
                          selected = ifelse(input$use_add == 1, "add compound", "parameters")
        )
    })
    
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
                                  "description"=parameter_names[input$cparams_select], 
                                  "value"=input$cparams_value, 
                                  "MC 2.5%" = NA, 
                                  "MC mean" = NA, 
                                  "MC 97.5%" = NA)
            # browser
            if(exists("custom_param_values")) {
                #remove the last existing value if it was in there
                custom_param_values <<- custom_param_values[custom_param_values$parameter != input$cparams_select,]
                custom_param_values <<- rbind(custom_param_values, 
                                              newdata
                )
            } else {
                custom_param_values <<- newdata
            }
        
        #clean the inputs
        updateNumericInput(session, "cparams_value", value = 0)
    })
    # custom_param_values <- data.frame("parameter"=c(), "description"=c(), "value"=c(), 
                                      # "MC 2.5%"=c(), "MC mean"=c(), "MC 97.5%"=c())
    output$custom_param_table <- renderTable({
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
        # browser()
        }
        inits <- do.call(parameterize_pbtk, param_list)
        
        
        #update if user supplied custom values
        input$cparams_submit
        if(exists("custom_param_values") && nrow(custom_param_values) > 0) {
            torep <- custom_param_values$value
            names(torep) <- custom_param_values$parameter
            inits[names(torep)] <- torep
        }
        return(inits)
        
    })
        
    output$parameters_df <- renderTable({
        # 
        ww <- parameters()
        parameter_df <- data.frame("parameter"=names(ww), "description" = parameter_names[names(ww)], "value"=unlist(ww))
        
        # if(input$output_type == "single")
        if(input$output_type == "mc") {
            #expand the data frame with uncertainty info:
            df <- apply(do.call(rbind, lapply(results_mc()[["inits"]], unlist)), 2, function(x) {
                c("mean"=mean(x, na.rm=T), "lci"=quantile(x,.025, na.rm=T), "uci"=quantile(x,.975, na.rm=T))
            })
            parameter_df[["MC 2.5%"]] <- df[2,]
            parameter_df[["MC mean"]] <- df[1,]
            parameter_df[["MC 97.5%"]] <- df[3,]
            
        }
            return(parameter_df)
            
    })
    
    results <- reactive({
        solve_list <- list(
            parameters = parameters(),
            chem.name=input$compound,
            plots=FALSE,
            suppress.messages = TRUE,
            output.units = input$solve.output.units,
            iv.dose = input$solve.iv.dose,
            tsteps = input$solve.tsteps,
            days = input$solve.days
        )
        if(input$dose_type == "daily dose")
            solve_list$daily.dose=input$solve.daily.dose
        if(input$dose_type == "per dose + doses/day") {
            solve_list$doses.per.day=input$solve.doses
            solve_list$dose=input$solve.dose
        }
        
        if(input$output_type == "single") 
            return(do.call(solve_pbtk, solve_list))
        

    })
    
    results_mc <- eventReactive(input$run, {
        if(input$output_type == "mc") {
            pbtk_mclist <- list("inits"=list(), "pbtk_result"=list(), "halflife"=list(), "AUC"=list(), "Cmax"=list())
            param_list <- list("chem.cas"=NULL,"chem.name" = input$compound, 
                               "species" = input$species, "default.to.human" = F,
                               "tissuelist" = list(liver=c("liver"), 
                                                   kidney=c("kidney"), 
                                                   lung=c("lung"), 
                                                   gut=c("gut")),
                               "force.human.clint.fub" = F, "clint.pvalue.threshold" = 0.05, 
                               "monte.carlo" = TRUE, 
                               "monte.carlo.cv" = mc_cv(),
                               "monte.carlo.log" = input$mc_use_log)
            if(input$use_cas) {
                param_list$chem.cas <- input$cas
                param_list$chem.name <- NULL
            }
            solve_list <- list(
                chem.name=input$compound,
                plots=FALSE,
                suppress.messages = TRUE,
                output.units = input$solve.output.units,
                iv.dose = input$solve.iv.dose,
                tsteps = input$solve.tsteps,
                days = input$solve.days
            )
            if(input$dose_type == "daily dose")
                solve_list$daily.dose=input$solve.daily.dose
            if(input$dose_type == "per dose + doses/day") {
                solve_list$doses.per.day=input$solve.doses
                solve_list$dose=input$solve.dose
            }
            
            withProgress(message = "Calculation in progress", 
                         detail = paste0("Performing ", input$nSimulations, " Monte Carlo estimates"),
                         for(i in 1:input$nSimulations) {
                            incProgress(1/input$nSimulations)
                            inits <- do.call(parameterize_pbtk, param_list)
                            solve_list$parameters <- inits
                            pbtk_mclist[["inits"]][[i]] <- inits
                            sol <- do.call(solve_pbtk, solve_list)
                            
                            pbtk_mclist[["pbtk_result"]][[i]] <- sol
                            times <- sol[,"time"]
                            x <- sol[,"Cplasma"]
                            wmax <- which.max(x)
                            wmin <- which(x < (max(x)/2))
                            pbtk_mclist[["halflife"]][[i]] <- times[min(wmin[wmin > wmax])] - times[wmax]
                            pbtk_mclist[["AUC"]][[i]] <- llTrapAUC(times, x)
                            pbtk_mclist[["Cmax"]][[i]] <- max(x)
                            
                         }, 
                         min = 0, max = 1)
            
            return(pbtk_mclist)
    }      
    })
    
    endpoints <- reactive({
        # ww <- c("Cplasma", paste0("C", input$compartments), "Crest", "Ametabolized", "Atubules", "Agutlumen")
        ww <- c("Cplasma", paste0("C", c("lung", "kidney", "gut", "liver")), "Crest", "Ametabolized", "Atubules", "Agutlumen")
        # names(ww) <- c("Plasma", input$compartments, "rest", "metabolized", "tubules", "gut lumen")
        names(ww) <- c("Plasma", c("lung", "kidney", "gut", "liver"), "rest", "metabolized", "tubules", "gut lumen")
        return(ww)
    })
    results_mc_df <- reactive({
        lci <- (1-input$display_ci)/2
        uci <- 1 - (1-input$display_ci)/2
        res <- apply(abind::abind(results_mc()[["pbtk_result"]], along=3), c(1,2), function(x) {
            c("mean"=mean(x, na.rm=T), "lci"=quantile(x,lci, na.rm=T), "uci"=quantile(x,uci, na.rm=T))
        })
        select <- c("time", endpoints())
        if(!is.null(select)) res <- res[,,select] #you can select some of the columns only
        return(res)
    })
    
    output$results_plot <- renderPlot({
        # 
        res <- results_mc_df()
        if(input$output_type == "mc") {
            timevar <- res["mean",,"time"]
            par(mfrow=c(ceiling(dim(res)[3]/3),3))
            for(cd in 1:dim(res)[3]) {
                plot(res["mean",,cd] ~ timevar, type="l", 
                     xlab="time (days)", ylab=dimnames(res)[3][[1]][cd])
                polygon(c(timevar, rev(timevar)), 
                        c(res[2,,cd], rev(res[3,,cd])), 
                        col="gray", border=NA)
                lines(res["mean",,cd] ~ timevar, type="l", lwd=1.2)
            }
            
        }
        if(input$output_type == "single")  {
            par(mar=c(1,1,1,1))
            plot(results())
        }
    })
    
    output$choose_plot_ui <- renderUI({
        selectInput("choose_plot", "Choose plot for detailed display", endpoints())
    })
    
    output$results_plot_single <- renderPlot({
        
        if(!is.null(input$choose_plot)) {
            if(input$output_type == "mc") {
                res <- results_mc_df()
                timevar <- res["mean",,"time"]
                cd <- which(dimnames(res)[[3]] == input$choose_plot)
                # browser()
                plot(res["mean",,cd] ~ timevar, type="l", xlab="time (days)", ylab=dimnames(res)[3][[1]][cd])
                polygon(c(timevar, rev(timevar)), c(res[2,,cd], rev(res[3,,cd])), col="gray", border=NA)
                lines(res["mean",,cd] ~ timevar, type="l", lwd=1.2)
            }
            if(input$output_type == "single") {
                res <- results()
                cd <- which(colnames(res) == input$choose_plot)
                plot(res[,cd] ~ res[,"time"], type="l", xlab="time (days)", ylab=endpoints()[cd])
            }
        }
    })
    
    output$results_plot_ui <- renderUI({
        if(input$output_type == "mc")
            plotOutput("results_plot")
            # plotOutput("results_plot", height=600, width= 600)
        if(input$output_type == "single")
            plotOutput("results_plot", height=200, width= 600)
        
    })
    
    output$results_numerical <- renderTable({
        if(input$output_type == "mc") {
            lci <- (1-input$display_ci)/2
            uci <- 1 - (1-input$display_ci)/2
                
            df <- apply(data.frame(unlist(results_mc()[["halflife"]]), 
                                   unlist(results_mc()[["AUC"]]), 
                                   unlist(results_mc()[["Cmax"]])), 
                        2, function(x) {
                c("lci"=quantile(x,lci, na.rm=T), "mean"=mean(x, na.rm=T), "uci"=quantile(x,uci, na.rm=T)) 
            })    
                # browser()
            rownames(df) <- c(paste0(100*lci, "%"), "mean", paste0(100*uci, "%"))
            colnames(df) <- c("Half-life (Cplasma)", "AUC (Cplasma)", "Cmax (Cplasma)")
            # browser()
            return(df)
        }
        if(input$output_type == "single") {
            if(!is.null(results())) {
                times <- results()[,"time"]
                x <- results()[,"Cplasma"]
                wmax <- which.max(x)
                wmin <- which(x < (max(x)/2))
                tt <- times[min(wmin[wmin > wmax])] - times[wmax]
                return(data.frame("Cplasma half-life" = tt,
                                  "Cplasma Cmax" = max(x),
                                  "Cplasma AUC" = llTrapAUC(times, x)))
            }
        }
    }, rownames = TRUE, digits = 3)
    
    output$fileDownload <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".csv", sep="")
        },
        content = function(file) {
            # data <- ifelse(input$output_type == "single", results(), results_mc()[["pbtk_result"]][[1]])
            if(input$output_type == "single")
                data <- results()
            if(input$output_type == "mc")
                data <- results_mc_df()["mean",,]
            # browser()
            write.csv(data, 
                      file)
        },
        contentType='text/csv'
    )
})
