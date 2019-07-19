#' Single Cell Load Data UI
#'
#' @export
#' @return None
sc_loadDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar layout with input and output definitions ----
    sidebarPanel(
      radioButtons(
        ns("loadData"),
        label = "Data Type",
        c("Count Data Table" = 1, "10x Genomics Data" = 2)
      ),

      conditionalPanel(
        condition = "input.loadData == 1",
        ns = ns,



        radioButtons(
          ns("sep"),
          "Separator",
          choices = c(
            Comma = ",",
            Space = " ",
            Tab = "\t",
            Semicolumn = ";"
          )
        ),


        fileInput(
          ns("file1"),
          "Choose Counts File",
          multiple = FALSE,
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        )

      ),

      conditionalPanel(
        condition = "input.loadData == 2",
        ns = ns,
        shinyDirButton(
          ns('directoryButton'),
          'Select Directory',
          'Please select a folder'
        )

      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(verbatimTextOutput(ns("loadDataText"), placeholder = T))
  )
}


#' Single Cell Load Data Server
#'
#' Enables the upload of two types of data: 10x Genomics and Count data
#'
#' @export
#' @return counts - The count table to be used in the sc analysis
sc_loadData <- function(input, output, session) {
  counts <- reactiveValues()

  # Load 10x -----
  observeEvent(input$directoryButton, {
    volumes <- getVolumes()

    shinyDirChoose(input,
                   'directoryButton',
                   roots = volumes,
                   session = session)

    path1 <- reactive({
      return(print(parseDirPath(volumes, input$directoryButton)))
    })

    # show meta data table
    if (!is.null(path1)) {

        req(nchar(path1()) > 0)
        counts$countTable <- load10xData(path1(), session)

    }
  })

  # Load CountTable -----
  observeEvent(input$file1, {
    counts$countTable <- read.csv(input$file1$datapath,
                                  sep = input$sep,
                                  row.names = 1)

  })

  observe({
    output$loadDataText <- renderText({
      sprintf(
        " Number of cells is %s;\n Number of genes is %s;",
        ncol(counts$countTable),
        nrow(counts$countTable)
      )
    })

  })

  return(counts)
}


#' Load 10x Data with TryCatch
#'
#' Simply adds TryCatch statement to the function supplied by Seurat
#'
#'
#' @param path The directory containing the 10X genomics data
#'
#' @export
#' @return
load10xData <- function(path, session) {
  counts <- tryCatch(
    {
      counts <- Read10X(data.dir = path)
    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data format error",
        text = "Ensure that a folder containing 10X Genomics data was appropriately chosen",
        type = "error"
      )
      return()
    }
  )
  return(counts)
}
