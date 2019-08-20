#' Single Cell Compare Data UI
#'
#' @export
#' @return None
sc_compUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarPanel(

      h4("Generate Venn Diagram"),

      numericInput(
        ns("compMpctInput"),
        label = "Minimum Cell Fraction of Genes",
        value = 0.01,
        min = 0,
        max = 0.5
      ),

      numericInput(
        ns("compFCInput"),
        label = "Fold-Change Threshold >",
        value = 2,
        min = 1,
        max = 10
      ),
      numericInput(
        ns("compPvalueInput"),
        label = "Adjusted P-value Threshold <",
        value = 0.05,
        min = 0.000001,
        max = 0.5
      ),
      actionButton(ns("comparisonButton"), "Run DE Pipelines"),

      tags$hr(),

      conditionalPanel(
        condition = "input.comparisonButton > 0",
        ns = ns,

        selectInput(
          ns("selectIntersect"),
          label = ("Select intersect of interest"),
          choices = list(
            "Wilcoxon (only)" = "a1",
            "Wilcoxon & t-test" = "a2",
            "t-test (only)" = "a3",
            "DESeq2 x MAST" = "a4",
            "All Three Methods" = "a5",
            "t-test & MAST" = "a6",
            "MAST (only)" = "a7"
          ),
          selected = NULL
        ),

        downloadButton(ns("downloadIntersect"), "Download Intersect")

      )
    ),

    mainPanel(
      htmlOutput(ns("helpCompInfo")),
      plotOutput(ns("comparsionPlot"), width = "800px", height = "500px")
      )
  )
}



#' SC Compare Data Server
#'
#' @param finData Clustered data results
#' @export
#' @return None
sc_comp <- function(input, output, session, finData) {
  comp <- reactiveValues()

  output$helpCompInfo <- renderUI({
    if(input$comparisonButton == 0){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px;'>
        <p style ='text-align: center'><b>
        This tab supplies users with an option to assess
        the agreement between the different DE analysis packages.</b> </p> <br>

        Prior to running the pipeline,
        users can pre-filter genes according to: <br>
        Fold-change, adj. P-value threshold,
        and genes expressed in a minimum fraction of cells. <br> <br>

        Once the pipeline is finished a Venn Diagram
        with the overlap between selected DE methods is returned.
        Each overlap(intersect) can then be downloaded <br> <br>

        <i>Note that the procedure runs 4 subsequent DE analysis pipelines,
           as such it is rather time-consuming.</i> </div>" )
    } else {
      HTML("")
    }
  })

  observeEvent(input$comparisonButton, {
    #* quite different

    show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))

    comp$xlist <-
      sc_getAllDE(
        finData$finalData,
        input$compMpctInput,
        input$compFCInput,
        input$compPvalueInput
      )

    grid.newpage()

    comp$plot <- (plotAllVenn(comp$xlist))

    output$comparsionPlot <- renderPlot({
      grid.draw(comp$plot)
    })

    hide_waiter()

  })

  observeEvent(input$selectIntersect, {
    if (!is.null(comp$xlist)) {
      comp$intersect <- getIntersect(comp$xlist, input$selectIntersect)
    }

  })

  output$downloadIntersect <- downloadHandler(
    filename = function() {
      paste("output/PackageComparison_intersect",
            input$selectIntersect,
            ".csv",
            sep = "")
    },
    content = function(file) {
      data <- comp$intersect

      write.csv(data, file, row.names = FALSE)
    }
  )
}



#' SC Generate Data required to compare DE Method Results
#'
#' @param data Clustering results
#' @param mPCT Test Genes detected in a minimum fraction of min.pct cells
#' @param fc Fold-change threshold
#' @param pValue p-value threshold
#' @export
#' @return Returns a list with DE genes according to the different methods
sc_getAllDE <- function(data, mPCT, fc, pValue) {
  fc <- log(fc)

  x1 <-
    FindAllMarkers(
      data,
      test.use = "wilcox",
      min.pct = mPCT,
      logfc.threshold = fc
    )
  x2 <-
    FindAllMarkers(
      data,
      test.use = "t",
      min.pct = mPCT,
      logfc.threshold = fc
    )

  x3 <-
    FindAllMarkers(
      data,
      test.use = "MAST",
      min.pct = mPCT,
      logfc.threshold = fc
    )


  x1_sig <- subset(x1, p_val_adj < pValue)
  x2_sig <- subset(x2, p_val_adj < pValue)
  x3_sig <- subset(x3, p_val_adj < pValue)


  list <-
    list(rownames(x1_sig),
         rownames(x2_sig),
         rownames(x3_sig))

  names(list) <-
    c("Wilcoxon", "T-test", "MAST")

  return(list)
}
