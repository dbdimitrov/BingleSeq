#' Initiate GUI
#'
#'
#' @export
#' @return None
ui <-  tagList(
  useShinyjs(),
  use_waiter(),
      navbarPage(id = "mainPage",
                 title = "BingleSeq",

                 tabPanel(title = "Choose App",
                          value = "startApp",

                          radioGroupButtons("chooseApp",
                                            h4("Pick Data Type"),
                                            size = 'lg',
                                            choices =  c("scRNA-Seq Analysis" = 1,
                                                         "Bulk RNA-Seq Analysis" = 2)),

                          actionBttn(("launch_app"), "Launch App", size = "md", style = "material-flat")

                 )
    )
)
