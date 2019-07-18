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

                      radioButtons("chooseApp",
                                   "Pick Data Type",
                                   choices =  c("Single Cell" = 1, "Bulk" = 2)),

                      actionButton(("launch_app"), "Launch App"))

  )
)
