
library(shiny)
library(httk)

nSimulations <- 1000
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
    "Vvenc" = "Volume of the veins per kg body weight, L/kg BW."
)

shinyServer(function(input, output) {

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
    
    inits_mean <- parameterize_pbtk(chem.cas=NULL,chem.name = "Bisphenol-A", species = "Human", default.to.human = F,
                                    tissuelist = list(liver=c("liver"), kidney=c("kidney"), lung=c("lung"), gut=c("gut")),
                                    force.human.clint.fub = F, clint.pvalue.threshold = 0.05, monte.carlo=FALSE)
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
        do.call(parameterize_pbtk, param_list)
    })
        
    output$parameters_df <- renderTable({
        ww <- parameters()
        parameter_df <- data.frame("parameter"=names(ww), "description" = parameter_names, "value"=unlist(ww))
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
            pbtk_mclist <- list("inits"=list(), "pbtk_result"=list(), "halflife"=list())
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
                         }, 
                         min = 0, max = 1)
            
            return(pbtk_mclist)
    }      
    })
    
    endpoints <- reactive({
        ww <- c("Cplasma", paste0("C", input$compartments), "Crest", "Ametabolized", "Atubules", "Agutlumen")
        names(ww) <- c("Plasma", input$compartments, "rest", "metabolized", "tubules", "gut lumen")
        return(ww)
    })
    results_mc_df <- reactive({
        res <- apply(abind::abind(results_mc()[["pbtk_result"]], along=3), c(1,2), function(x) {
            c("mean"=mean(x, na.rm=T), "lci"=quantile(x,.025, na.rm=T), "uci"=quantile(x,.975, na.rm=T))
        })
        # timevar <- res["mean",,"time"]
        # select <- dimnames(res)[3][[1]][-1]
        select <- c("time", endpoints())
        if(!is.null(select)) res <- res[,,select] #you can select some of the columns only
        return(res)
    })
    
    output$results_plot <- renderPlot({
        if(input$output_type == "mc") {
            res <- results_mc_df()
            timevar <- res["mean",,"time"]
            par(mfrow=c(ceiling(dim(res)[3]/3),3))
            for(cd in 1:dim(res)[3]) {
                plot(res["mean",,cd] ~ timevar, type="l", xlab="time", ylab=dimnames(res)[3][[1]][cd])
                polygon(c(timevar, rev(timevar)), c(res["lci.2.5%",,cd], rev(res["uci.97.5%",,cd])), col="gray", border=NA)
                lines(res["mean",,cd] ~ timevar, type="l", lwd=1.2)
            }
        }
        if(input$output_type == "single") 
            plot(results())
    })
    
    output$choose_plot_ui <- renderUI({
        selectInput("choose_plot", "Choose plot for detailed display", endpoints())
    })
    
    output$results_plot_single <- renderPlot({
        if(input$output_type == "mc") {
            res <- results_mc_df()
            timevar <- res["mean",,"time"]
            cd <- which(endpoints() == input$choose_plot)
            plot(res["mean",,cd] ~ timevar, type="l", xlab="time", ylab=dimnames(res)[3][[1]][cd])
            polygon(c(timevar, rev(timevar)), c(res["lci.2.5%",,cd], rev(res["uci.97.5%",,cd])), col="gray", border=NA)
            lines(res["mean",,cd] ~ timevar, type="l", lwd=1.2)
        }
        if(input$output_type == "single") {
            res <- results()
            cd <- which(endpoints() == input$choose_plot)
            # browser()
            
            plot(res[,cd] ~ res[,"time"], type="l", xlab="time", ylab=endpoints()[cd])
            
        }
    })
    
    

})
