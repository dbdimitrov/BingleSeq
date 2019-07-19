#' Initiate GUI
#'
#'
#' @export
#' @return None
ui <-  tagList(
  useShinyjs(),
  navbarPage(id = "mainPage",
             title = "BingleSEQ",

             tabPanel(title = "Choose App",
                      value = "startApp",

                      h5("Choose directory to save additional files"),

                      shinyDirButton('saveFilesDirButton', 'Select Directory', 'Please select a folder to save output'),


                      radioButtons("chooseApp",
                                   "Pick Data Type",
                                   choices =  c("Single Cell" = 1, "Bulk" = 2)),


                      actionButton(("launch_app"), "Launch App"))

  )
)
