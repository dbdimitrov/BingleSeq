#' Single Cell Normalize Tab UI
#'
#' @export
#' @return None
sc_normUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(
      h4("Normalize the Data"),
      selectInput(
        ns("normalizeCombo"),
        label = "Select Normalization method",
        choices = list(
          "Log Normalize" = "LogNormalize",
          "Relative counts" = "RC"
        ),
        selected = 1
      ),

      numericInput(
        ns("scaleNumInput"),
        label = "Set Scale Factor for normalization",
        min = 1000,
        value = 10000
      ),

      tags$hr(),

      h4("Find Variable Features"),

      selectInput(
        ns("varianceCombo"),
        label = "Variance Estimation Method",
        choices = list(
          "VST" = "vst",
          "Mean Variance Plot" = "mvp",
          "Dispersion" = "disp"
        ),
        selected = 1
      ),

      numericInput(
        ns("varianceFeatInput"),
        label = "Set FeatureNo for Variance Estimation",
        min = 100,
        value = 2000
      ),

      tags$hr(),

      actionButton(ns("normalizeButton"), label = "Find Variable Genes")
    ),

    # Main panel for displaying outputs ----
    mainPanel(
      htmlOutput(ns("helpNormInfo")),
      plotlyOutput(ns("normalizePlot"), width = "1280px", height = "720px"),
      textOutput(ns("violionPlotInfo"))
      )
  )
}


#' Single Cell Normalize Tab Server
#'
#' @param filtData Reactive value containing the suerat Object with filtered data
#'
#' @export
#' @return Returns a reactive value containing Seurat object with normalized data
sc_norm <- function(input, output, session, filtData) {
  norm <- reactiveValues()

  output$helpNormInfo <- renderUI({
    if(is.null(norm$normalizedData)){
    HTML("<div style='border:2px solid blue; padding-top: 8px; padding-bot: 8px; font-size: 14px;
      border-radius: 10px;'>
    <p style='text-align: center'><b> This tab enables the normalization of the data. </b> </p> <br>

    For integer counts use the 'Log Normalize' method,
    while for relative counts (e.g. FPMKs, CPM) use 'Relative Counts' normalization. <br>
    10,000 Scale Factor should be appropriate for almost any type of data,
    besides when working with CPM as 1,000,000 should be used in this case. <br> <br>

    Simultaneously with Normalization, the Most Variable Features in the dataset will be identified.
    The number of Variable Features is set with 'FeatureNo'
    and variance is estimated with a method of choice. <br>
    Following this process a Variance plot with the Most
    Variable Features is displayed. <br> <br>

    <i>Note that: the identified Most Variable Features
         are relavant only for unsupervised clustering with Seurat.</i> </div>")
    } else{
      HTML("")
    }

  })

  ### Normalization ------
  observeEvent(input$normalizeButton, {

    waiter_show(html = tagList(spin_folding_cube(), h2("Loading ...")))
    norm$normalizedData <-
      normalizeSeurat(
        filtData$filteredData,
        input$normalizeCombo,
        input$scaleNumInput,
        input$varianceCombo,
        input$varianceFeatInput
      )


    if (!is.null(norm$normalizedData)) {
      norm$variancePlot <- VariableFeaturePlot(norm$normalizedData)


    }

    output$normalizePlot <- renderPlotly({
      ggplotly(norm$variancePlot)
      style(norm$variancePlot, text = VariableFeatures(norm$normalizedData))
    })

    output$violionPlotInfo <- renderText({
      "Note: Highly Variable Features are shown in red (Hover for Gene names)"
    })

    waiter_hide()

  })

  return(norm)
}


#' Single Cell Normalize function
#'
#' @param s_object Suerat Object with filtered data
#' @param normalizeMet The normalization method to be used
#' @param scaleF The Scale Factor
#' @param varianceMet The variance estimation method to be used
#' @param nfeat The number of features to be used in estimating variance
#' @export
#' @return Returns a Seurat object with normalized data
normalizeSeurat <-
  function(s_object,
           normalizeMet,
           scaleF,
           varianceMet,
           nfeat) {
    normalized_object <-
      NormalizeData(s_object,
                    normalization.method = normalizeMet,
                    scale.factor = scaleF)

    if (startsWith(varianceMet, "mvp")) {
      normalized_object <-
        FindVariableFeatures(normalized_object,
                             selection.method = varianceMet)
    } else{
      normalized_object <-
        FindVariableFeatures(normalized_object,
                             selection.method = varianceMet,
                             nfeatures = nfeat)
    }

    return(normalized_object)
  }
