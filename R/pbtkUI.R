#this launches a bundled-in Shiny application in a separate window
pbtkUI <- function() {
    shiny::runApp(system.file('shiny', package="httk"), display.mode = "normal")
}