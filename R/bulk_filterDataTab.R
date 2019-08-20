#' Bulk Filter Tab UI
#'
#' @export
#' @return None
bulk_filterDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarPanel(
      h4("Filter Counts Table"),

      # Select Filter Method
      selectInput(
        ns("selectFilter"),
        label = "Select Filter Method",
        choices = list(
          "CPM" = 1,
          "Median" = 2,
          "Max" = 3
        ),
        selected = 1
      ),

      tags$hr(),

      numericInput(ns("filterValue"), label = "Filter Genes with CPM below", value = 1),

      conditionalPanel(
        condition = "input.selectFilter == 1",
        ns = ns,
        numericInput(ns("sampleNoInput"), label = "in at least # samples", value = 3)

      ),

      hr(),

      fluidRow(column(3, verbatimTextOutput(
        ns("filterValue")
      ))),

      actionButton(ns("filterButton"), label = "Filter Data")

    ),

    # Main panel for displaying outputs ----
    mainPanel(fluidRow(
      column(
        width = 4,
        offset = 0,
        tags$h4("Raw Data Summary"),
        DT::dataTableOutput(ns("prefilterTable"))

      ),

      column(
        width = 4,
        offset = 0,
        tags$h4("Filtered Data Summary"),
        DT::dataTableOutput(ns("postfilterTable"))

      )
    ),
    fluidRow(
      column(width = 4,
             offset = 0,
             plotOutput(ns("preFiltHist"))),

      column(width = 4,
             offset = 0,
             plotOutput(ns("postFiltHist")))
    ))
  )
}


#' Bulk Filter Tab Server
#'
#' @param counts The loaded Count Table
#'
#' @export
#' @return Returns a reactive value with the filtered Count Table
bulk_filterData <- function(input, output, session, counts) {
  filt <- reactiveValues()

  observeEvent(input$filterButton, {
    filt$filteredCounts <-
      filterFunction(
        counts$countTable,
        as.numeric(input$selectFilter),
        as.numeric(input$filterValue),
        as.numeric(input$sampleNoInput)
      )

    output$postfilterTable <- DT::renderDataTable(DT::datatable(
      generateSummary(filt$filteredCounts, session),
      options = list(paging = FALSE, searching = FALSE),
      rownames = FALSE
    ))

    output$postFiltHist <- renderPlot({
      qcHist(filt$filteredCounts)
    })

    # Used to generate DE Tab only when generateSummary is OK
    filt$correctFormat <- TRUE
  })


  observeEvent (input$selectFilter , {
    if (input$selectFilter == 1) {
      updateNumericInput(session,
                         "filterValue",
                         label = "Filter Genes with CPM below",
                         value = 1)
    } else if (input$selectFilter == 2) {
      updateNumericInput(session,
                         "filterValue",
                         label = "Filter Genes with Row Median below",
                         value = 10)
    } else{
      updateNumericInput(session,
                         "filterValue",
                         label = "Filter Genes with Row Maximum below",
                         value = 10)
    }
  })

  return(filt)
}


#' Filter Function
#'
#' @param countTable The loaded Count Table
#' @param method Filter Option (CPM, MAX, Median)
#' @param value Keep only genes with counts > than this value
#'
#' @export
#' @return Returns the filtered Count Table
filterFunction <- function(data, method, value, sampleNo) {
  countTable  <- data[, -1]
  rownames(countTable) <- data[, 1]


  #filter by CPM
  if (method == 1) {
    keep <- rowSums(edgeR::cpm.default(countTable) > value) >= sampleNo
    countTable <- countTable[keep, ]

    #filter by Median
  } else if (method == 2) {
    countTable <-
      subset(countTable, apply(countTable, 1, median) >= value)

    #filter by Max
  } else if (method == 3) {
    countTable <- subset(countTable, apply(countTable, 1, max) >= value)

  }


  IDs <- rownames(countTable)

  rownames(countTable) <- NULL

  countTable <- cbind(IDs, countTable)

  return(countTable)
}


#' QC Histogram Function
#'
#' @param data Filtered/Unfiltered Count Table
#'
#' @export
#' @return returns a QC Histogram
qcHist <- function(data) {

  log10_data <- log10(data[, 2:ncol(data)])

  p <-
    hist(
      rowSums(log10_data[, 2:ncol(log10_data)]),
      col = "grey",
      main = "",
      xlab = "log10 Total Counts per Gene",
      ylab = "Frequency",
      breaks = 50,
      xlim = range(0:max(rowSums(log10_data)))
    )

  return(p)
}
