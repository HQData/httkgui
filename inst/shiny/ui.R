
# compartment_names <- c("lung", "kidney", "gut", "liver")
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

shiny::shinyUI(fluidPage(

  # Application title
  titlePanel("httk UI for PBTK models"),

  sidebarLayout(
    sidebarPanel(
        fluidRow(
            column(7, checkboxInput("use_cas", "Use CAS instead of compound name")),
            column(5, checkboxInput("use_add", "Add a new compound"))
        ),
      conditionalPanel("input.use_add == 0",
          conditionalPanel("input.use_cas == 0", selectizeInput("compound", "Compound name", 
                                                                httk::chem.physical_and_invitro.data$Compound)),
          conditionalPanel("input.use_cas == 1", selectizeInput("cas", "CAS", 
                                                                httk::chem.physical_and_invitro.data$CAS))
      ),
      conditionalPanel("input.use_add == 1", "Please define the compound to use in separate tab on the right"),
      selectInput("species", "Species", c("Human", "Rat", "???")),
      # checkboxGroupInput("compartments", "Compartments of interest", compartment_names, selected=compartment_names),
      # specify options for solve_pbtk: daily.dose, dose, doses.per.day, species, iv.dose, output.units, tsteps, days
      radioButtons("dose_type", "Dose specification", c("daily dose", "per dose + doses/day"), inline=TRUE),
      conditionalPanel('input.dose_type == "daily dose"', 
        numericInput("solve.daily.dose", "Single dose (mg/kg BW)", 1)
      ),
      conditionalPanel('input.dose_type == "per dose + doses/day"', 
        numericInput("solve.doses.per.day", "Doses per day", 0),
        numericInput("solve.dose", "Per dose/day (mg/kg BW)", 0)
      ),
      checkboxInput("solve.iv.dose", "IV dose (default = oral)", 0),
      radioButtons("output_type", "Single simulation or Monte Carlo?", c("single"="single", "Monte Carlo"="mc"), inline=TRUE),
      selectInput("solve.output.units", "Output units", c("mg/L", "mg", "umol", "uM"), selected="uM"),
      fluidRow(
        column(6, numericInput("solve.tsteps", "time steps / hour", 4)),
        column(6, numericInput("solve.days", "Simulation length (days)", 10))
        ),
      conditionalPanel("input.output_type == 'mc'", actionButton("run", "Solve PBPK model"))
    ),

    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(id="main_panel",
      tabPanel("parameters", 
        h3("PBTK model parameter values"),
        checkboxInput("custom_params", "Check here to manually change parameter values", 0, width=500),
        conditionalPanel("input.custom_params == 1",
                         em("To erase user-defined values, please restart the application"),
                         
                         fluidRow(
                             column(4, selectInput("cparams_select", "Parameter name", names(parameter_names))),
                             column(2, numericInput("cparams_value", "Mean value", 0)),
                             
                             conditionalPanel("input.output_type == 'mc'",
                                              column(2, numericInput("cparams_lci", "2.5% value", 0)),
                                              column(2, numericInput("cparams_uci", "97.5% value", 0))
                             )
                         ),
                         actionButton("cparams_submit", "Submit values"),
                         tableOutput("custom_param_table"),
                         h3("Original parameter values"),
                         em("These values will be replaced by user-defined values.")
        ),
        
        tableOutput("parameters_df")
      ),
      tabPanel("add compound",
               conditionalPanel("input.use_add == 1", 
                                h3("Define new compound and press submit when all inputs ready"),
                                column(4,
                                    textInput("add_compound", "Compound name", ""),
                                    textInput("add_cas", "Compound CAS", ""),
                                    textInput("add_reference", "Reference(s)", value=""),
                                    HTML("Species can be selected in the input panel on the left <br><br>"),
                                    HTML("<em>If data is not available, please check the NA box instead of inputting a value.</em> <br><br>"),
                                    actionButton("add_submit", "Submit")
                                ),
                                column(4, 
                                    numericInput("add_mw", "MW", value=0),
                                    checkboxInput("add_mw_na", "MW NA", value=0),
                                    numericInput("add_logp", "LogP", value=0),
                                    checkboxInput("add_logp_na", "LogP NA", value=0),
                                    numericInput("add_funbound", "Funbound plasma", value=0),
                                    checkboxInput("add_funbound_na", "Funbound plasma NA", value=0),
                                    numericInput("add_fgutabs", "Fgutabs", value=0),
                                    checkboxInput("add_fgutabs_na", "Fgutabs NA", value=0),
                                    textInput("add_pka_donor", "pKa donor", value=""),
                                    checkboxInput("add_pka_donor_na", "pKa donor NA", value=0),
                                    textInput("add_pka_accept", "pKa accept", value=""),
                                    checkboxInput("add_pka_accept_na", "pKa accept NA", value=0)
                                ),
                                column(4, 
                                    numericInput("add_clint", "Clint", value=0),
                                    checkboxInput("add_clint_na", "Clint NA", value=0),
                                    numericInput("add_kts", "KTS", value=0),
                                    checkboxInput("add_kts_na", "KTS NA", value=0),
                                    numericInput("add_fr", "FR", value=0),
                                    checkboxInput("add_fr_na", "FR NA", value=0),
                                    numericInput("add_vmax", "Vmax", value=0),
                                    checkboxInput("add_vmax_na", "Vmax NA", value=0),
                                    numericInput("add_km", "Km", value=0),
                                    checkboxInput("add_km_na", "Km NA", value=0)
                                )
                                
                                
                                
               ),
               conditionalPanel("input.use_add == 0", "Please select 'add compound' on the left to manually enter compound data'")
      ),
      tabPanel("Monte Carlo settings", 
        conditionalPanel("input.output_type == 'mc'", 
                       h3("Define Monte Carlo simulation parameters"),
                       numericInput("nSimulations", "Number of draws", 50),
                       checkboxInput("mc_use_log", "Use log-normal distributions", 1),
                       h3("Coefficient of variation values"),
                       column(4,
                        sliderInput("cv.water",		'Total Body Water' 			, 0, 1, .3), 
                        sliderInput("cv.plasma",	    'Plasma Volume' 			, 0, 1, .3),
                        sliderInput("cv.cardiac",	'Cardiac Output' 			, 0, 1, .3),
                        sliderInput("cv.bw",		    'Average BW' 				, 0, 1, .16)
                       ), column(4,
                       sliderInput("cv.tpp", 		'Total Plasma Protein' 		, 0, 1, .14),
                       sliderInput("cv.albumin",	'Plasma albumin' 			, 0, 1, .1),
                       sliderInput("cv.a1agp",		'Plasma a-1-AGP' 			, 0, 1, .3),
                       sliderInput("cv.hematocrit", 'Hematocrit'				, 0, 1, .3)
                       ), column(4,
                       sliderInput("cv.urine",		'Urine' 					, 0, 1, .3),
                       sliderInput("cv.bile",		'Bile' 						, 0, 1, .3),
                       sliderInput("cv.gfr",		'GFR' 						, 0, 1, .3),
                       sliderInput("cv.abt",		'Average Body Temperature' 	, 0, 1, 0)
                       )
          ),
          conditionalPanel("input.output_type == 'single'",
                           "Please select Monte Carlo mode on the left to define parameter variability.")
      ),
      tabPanel("results", 
        h3("Presentation of results of PBTK model"),
        conditionalPanel("input.output_type == 'mc'", 
            sliderInput("display_ci", "Uncertainty interval", min = 0, max = 1, value = 0.95, step = .05)),
        htmlOutput("choose_plot_ui"),
        plotOutput("results_plot_single"),
        fluidRow(
            column(6, 
                h3("Cmax, AUC, Half-life calculations"),
                tableOutput("results_numerical")
                # plotOutput("results_plot")
            ),
            column(6,
                h3("Downloading model results"),
                conditionalPanel("input.output_type == 'mc'", 
                                 "Values for Monte Carlo provided as mean over all simulations"),
                                 # selectInput("download_choice_mc", 
                                             # "Choose dataset to download", 
                                             # c("mean estimates with 95% intervals"))),
                # actionButton("download", "Download dataset as .csv file")
                downloadLink("fileDownload", "Download model solution (.csv file)")
            ))
        
        # uiOutput("results_plot_ui")
      )
    )
  ))
))
