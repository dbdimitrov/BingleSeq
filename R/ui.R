#' Initiate GUI
#'
#'
#' @export
#' @return None
ui <-  tagList(
  useShinyjs(),
  use_waiter(),
      navbarPage(id = "mainPage",
                 title = "BingleSEQ",

                 tabPanel(title = "Choose App",
                          value = "startApp",

                          radioGroupButtons("chooseApp",
                                            h4("Pick Data Type"),
                                            size = 'lg',
                                            choices =  c("scRNA-Seq Analysis" = 1,
                                                         "Bulk RNA-Seq Analysis" = 2)),


                          actionBttn(("launch_app"), "Launch App", size = "md", style = "material-flat")
                          # actionBttn(("launch_bulk"), "Bulk RNA-Seq Data Analysis", size = "md", style = "material-flat")

                          # h5("Choose directory to save additional files"),

                          # shinyDirButton('saveFilesDirButton', 'Select Directory', 'Please select a folder to save output'),

                 )
    )
)
