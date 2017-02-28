
library(shiny)

compartment_names <- c("lung", "kidney", "gut", "liver")


shinyUI(fluidPage(

  # Application title
  titlePanel("httk UI for PBTK models"),

  sidebarLayout(
    sidebarPanel(
      radioButtons("output_type", "Single simulation or Monte Carlo?", c("single"="single", "Monte Carlo"="mc"), inline=TRUE),
      checkboxInput("use_cas", "Use CAS instead of compound name"),
      conditionalPanel("input.use_cas == 0", selectizeInput("compound", "Compound name", httk::chem.physical_and_invitro.data$Compound)),
      conditionalPanel("input.use_cas == 1", selectizeInput("cas", "CAS", httk::chem.physical_and_invitro.data$CAS)),
      selectInput("species", "Species", c("Human", "Rat", "???")),
      checkboxGroupInput("compartments", "Compartments of interest", compartment_names, selected=compartment_names),
      # specify options for solve_pbtk: daily.dose, dose, doses.per.day, species, iv.dose, output.units, tsteps, days
          radioButtons("dose_type", "Dose specification", c("daily dose", "per dose + doses/day"), inline=TRUE),
          conditionalPanel('input.dose_type == "daily dose"', 
            numericInput("solve.daily.dose", "Daily dose (mg/kg BW)", 1)
          ),
          conditionalPanel('input.dose_type == "per dose + doses/day"', 
            numericInput("solve.dose", "Dose (mg/kg BW)", 0),
            numericInput("solve.doses.per.day", "Doses per day", 0)
          ),
          checkboxInput("solve.iv.dose", "IV dose", 0),
          selectInput("solve.output.units", "Output units", c("mg/L", "mg", "umol", "uM"), selected="uM"),
          numericInput("solve.tsteps", "time steps / hour", 4),
          numericInput("solve.days", "Simulation length (days)", 10),
      conditionalPanel("input.output_type == 'mc'", actionButton("run", "Solve PBPK model"))

    ),

    # Show a plot of the generated distribution
    mainPanel(tabsetPanel(
      tabPanel("parameters", 
        # HTML("Input parameters:"),
        tableOutput("parameters_df")
      ),
      tabPanel("Monte Carlo settings", conditionalPanel("input.output_type == 'mc'", 
                       numericInput("nSimulations", "Number of simulations", 100),
                       checkboxInput("mc_use_log", "Use log-normal distributions", 1),
                       h3("Coefficient of variation values"),
                       column(4,
                        sliderInput("cv.water",		'Total Body Water' 			, 0, 1, .3), 
                        sliderInput("cv.plasma",	    'Plasma Volume' 			, 0, 1, .3),
                        sliderInput("cv.cardiac",	'Cardiac Output' 			, 0, 1, .3),
                        sliderInput("cv.bw",		    'Average BW' 				, 0, 1, .16)
                       ),
                       column(4,
                       sliderInput("cv.tpp", 		'Total Plasma Protein' 		, 0, 1, .14),
                       sliderInput("cv.albumin",	'Plasma albumin' 			, 0, 1, .1),
                       sliderInput("cv.a1agp",		'Plasma a-1-AGP' 			, 0, 1, .3),
                       sliderInput("cv.hematocrit", 'Hematocrit'				, 0, 1, .3)
                       ),
                       column(4,
                       sliderInput("cv.urine",		'Urine' 					, 0, 1, .3),
                       sliderInput("cv.bile",		'Bile' 						, 0, 1, .3),
                       sliderInput("cv.gfr",		'GFR' 						, 0, 1, .3),
                       sliderInput("cv.abt",		'Average Body Temperature' 	, 0, 1, 0)
                       )
      )),
      tabPanel("results", 
        plotOutput("results_plot"),
        htmlOutput("choose_plot_ui"),
        plotOutput("results_plot_single")
      )
    )
  ))
))
