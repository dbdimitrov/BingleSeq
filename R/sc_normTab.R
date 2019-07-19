#' Single Cell Normalize Tab UI
#'
#' @export
#' @return None
sc_normUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(
      selectInput(
        ns("normalizeCombo"),
        label = "Select Normalization method",
        choices = list(
          "Log Normalize" = "LogNormalize",
          "Centered log ratio transformation" = "CLR",
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
      # 1e6 for rc

      tags$hr(),



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

      actionButton(ns("normalizeButton"), label = "Normalize Data")
    ),

    # Main panel for displaying outputs ----
    mainPanel(#verbatimTextOutput("normalizeText", placeholder = T),
      plotOutput(ns("normalizePlot")))
  )
}


#' Single Cell Normalize Tab Server
#'
#' @param filtData Reactive value containing the suerat Object with filtered data
#'
#' @export
#' @return Reactive value containing Seurat object with normalized data
sc_norm <- function(input, output, session, filtData) {
  norm <- reactiveValues()

  ### Normalization ------
  observeEvent(input$normalizeButton, {
    norm$normalizedData <-
      normalizeSeurat(
        filtData$filteredData,
        input$normalizeCombo,
        input$scaleNumInput,
        input$varianceCombo,
        input$varianceFeatInput
      )


    if (!is.null(norm$normalizedData)) {
      norm$variancePlot <- variancePlotSeurat(norm$normalizedData)

    }

    output$normalizePlot <- renderPlot({
      norm$variancePlot

    })

    # ggsave("figures/variancePlot.png", plot = norm$variancePlot, device = png(),
    #        width = 9, height = 6, limitsize = FALSE)

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
#' @return Seurat object with normalized data
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

    if (startsWith(varianceMet, "vst")) {
      normalized_object <-
        FindVariableFeatures(normalized_object,
                             selection.method = varianceMet,
                             nfeatures = nfeat)
    } else{
      normalized_object <-
        FindVariableFeatures(normalized_object, selection.method = varianceMet)
    }

    return(normalized_object)
  }

#' Single Cell Plot Variance Function
#'
#' @param s_object Suerat Object with filtered data
#' @export
#' @return Variance estimation plot with the ten most variable genes
variancePlotSeurat <- function(s_object) {
  # Identify the 10 most highly variable genes
  top10 <- head(VariableFeatures(s_object), 10)

  plot1 <- VariableFeaturePlot(s_object)
  plot2 <-
    LabelPoints(
      plot = plot1,
      points = top10,
      repel = T,
      xnudge = 0,
      ynudge = 0
    )

  return(plot2)
}
