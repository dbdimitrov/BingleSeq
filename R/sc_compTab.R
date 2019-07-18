#' Single Cell Compare Data UI
#'
#' @export
#' @return None
sc_compUI <- function(id) {
  ns <- NS(id)

  tagList(

    sidebarPanel(

      numericInput(ns("compMpctInput"), label = "Minimum Cell Fraction of Genes",
                   value = 0.25, min = 0, max = 0.5),

      numericInput(ns("compFCInput"), label = "Fold Change",
                   value = 2, min = 1, max = 10),
      numericInput(ns("compPvalueInput"), label = "P-value threshold",
                   value = 0.05, min = 0.000001, max = 0.5),
      actionButton(ns("comparisonButton"), "Generate Venn Diagram"),

      conditionalPanel(condition = "input.comparisonButton > 0",
                       ns = ns,

                       selectInput(ns("selectIntersect"), label =("Select intersect of interest"),
                                   choices = list("Negative Binomial (Only)" = "a1",
                                                  "Negative Binomial & MAST" = "a2",
                                                  "MAST (Only)" = "a3",
                                                  "Negative Binomial & Wilcoxon" = "a4",
                                                  "MAST & Negative Binomial & Wilcoxon" = "a5",
                                                  "All Four Statistics" = "a6",
                                                  "Negative Binomial & MAST & T-test" = "a7",
                                                  "MAST & T-test" = "a8",
                                                  "Wilcoxon (Only)" = "a9",
                                                  "Wilcoxon & MAST" = "a10",
                                                  "Wilcoxon & MAST & T-test" = "a11",
                                                  "Wilcoxon & Negative Binomial & T-test" = "a12",
                                                  "Negative Binomial & T-test" = "a13",
                                                  "T-test (Only)" = "a14",
                                                  "Wilcoxon & T-test" = "a15"),
                                   selected = NULL),

                       downloadButton(ns("downloadIntersect"), "Download Intersect")

      )
    ),

    mainPanel(
      plotOutput(ns("comparsionPlot"))
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

  observeEvent(input$comparisonButton, {   #* quite different

    comp$xlist <- sc_getAllDE(finData$finalData, input$compMpctInput,
                           input$compFCInput, input$compPvalueInput)

    grid.newpage()

    comp$plot <- (plotAllVenn(comp$xlist))

    output$comparsionPlot <- renderPlot({
      grid.draw(comp$plot)
    })

    # ggsave("figures/comparisonVenn.png", plot = comp$plot, device = png(),
    #        width = 12, height = 8, limitsize = FALSE)

  })

  observeEvent(input$selectIntersect, {

    if(!is.null(comp$xlist)){
      comp$intersect <- getIntersect(comp$xlist, input$selectIntersect)
    }

  })

  output$downloadIntersect <- downloadHandler(
    filename = function() {
      paste("output/PackageComparison_intersect",input$selectIntersect, ".csv", sep="")
    },
    content = function(file) {

      data <- comp$intersect;
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
#' @return A list with DE genes according to the different methods
sc_getAllDE <- function(data, mPCT, fc, pValue){

  fc <- log(fc)

  x1 <- FindAllMarkers(data, test.use = "wilcox", min.pct = mPCT, logfc.threshold = fc)
  x2 <- FindAllMarkers(data, test.use = "t", min.pct = mPCT, logfc.threshold = fc)
  data[["RNA"]]@counts<-as.matrix(data[["RNA"]]@counts)+1
  x3 <- FindAllMarkers(data, test.use = "negbinom", min.pct = mPCT, logfc.threshold = fc, slot = "counts")
  x4 <- FindAllMarkers(data, test.use = "MAST", min.pct = mPCT, logfc.threshold = fc)


  x1_sig <- subset(x1, p_val_adj < pValue)
  x2_sig <- subset(x2, p_val_adj < pValue)
  x3_sig <- subset(x3, p_val_adj < pValue)
  x4_sig <- subset(x4, p_val_adj < pValue)


  list <- list(rownames(x1_sig), rownames(x2_sig), rownames(x3_sig), rownames(x4_sig))

  names(list) <- c("Wilcoxon", "T-test", "Negative Binomial", "MAST")

  return(list)
}
