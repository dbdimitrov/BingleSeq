#' Single Cell Differential Expession Tab UI
#'
#' @export
#' @return None
sc_deUI <- function(id) {
  ns <- NS(id)
  tagList(# Sidebar panel for inputs ----
          sidebarPanel(tabsetPanel(
            id = ns("goSideTabSet"),
            tabPanel(
              title = "DGE test",

              selectInput(
                ns("dgeTestCombo"),
                label = "Select Test Type",
                choices = list(
                  "MAST" = "MAST",
                  "Wilcoxon Rank Sum test" = "wilcox",
                  "Student's T-test" = "t",
                  "ROC analysis" = "roc",
                  "Likelihood-ratio test(bimod)" = "bimod",
                  "Negative Binomial" = "Negative Binomial",
                  "Logistic Regression" = "LR"
                )
              ),

              numericInput(
                ns("logFCdgeInput"),
                label = "FC Threshold",
                min = 0,
                max = 10,
                value = 2
              ),

              numericInput(
                ns("pctdgeInput"),
                label = "Minimum Cell Fraction of Genes",
                min = 0,
                max = 0.5,
                value = 0.25
              ),

              actionButton(ns("dgeButton"), label = "Get Gene Markers"),


              conditionalPanel(
                condition = "input.dgeButton > 0",
                ns = ns,

                checkboxInput(ns("dgeClusterCheck"), h4("Show All Clusters"), TRUE),

                conditionalPanel(
                  condition = "!input.dgeClusterCheck",
                  ns = ns,

                  numericInput(
                    ns("dgeClustInput"),
                    label = "Cluster to Display",
                    min = 1,
                    max = 8,
                    value = 0
                  )
                )
              )
            ),

            tabPanel(
              title = "Plots",

              h4("Cluster Heatmap"),

              numericInput(
                ns("clustHeatInput"),
                label = "Genes to display",
                min = 1,
                value = 10
              ),

              actionButton(ns("dgeHeatButton"), label = "Generate Heatmap"),

              tags$hr(),

              h4("Choose Gene and Plot"),

              textInput(ns("geneNameInput"), "Enter Gene Name"),

              radioButtons(
                ns("dgePlotType"),
                label = "Plot Type",
                c(
                  "Violin Plot" = 1,
                  "Feature Plot" = 2,
                  "RidgePlot" = 3
                )
              ),

              actionButton(ns("dgePlotButton"), label = "Generate Plot")
            )
          )),

          # Main panel for displaying outputs ----
          mainPanel(tabsetPanel(
            id = ns("deMainTabSet"),
            tabPanel(title = "Table",

                     DT::dataTableOutput(ns("dgeTable"))),
            tabPanel(
              title = "Plot",
              value = "dePlotTab",

              plotOutput(ns("dgePlot")),
              downloadButton(ns("downloaddgePlot"), "Download Curret Plot")
            )
          )))
}

#' Single Cell Differential Expession Tab Server
#'
#' @param finData Reactive value containing a seurat object with clustered data
#'
#' @export
#' @return Diffenretial Expression data
sc_de <- function(input, output, session, finData) {
  de <- reactiveValues()

  ## Generate DE Data
  observeEvent(input$dgeButton, {
    # if(!is.null(finData$finalData)){

    de$markers <- FindAllMarkers(
      finData$finalData,
      test.use = input$dgeTestCombo,
      min.pct = input$pctdgeInput,
      logfc.threshold = log(input$logFCdgeInput)
    )



    # write.csv(de$markers, file="output/AllMarkerGenes.csv", row.names = FALSE)

    output$dgeTable <-
      DT::renderDataTable(if (input$dgeClusterCheck) {
        de$markers %>% datatable() %>%
          formatSignif(columns = c(1:2, 5), digits = 4)
      } else{
        de$markers[de$markers$cluster == input$dgeClustInput,] %>% datatable() %>%
          formatSignif(columns = c(1:2, 5), digits = 4)
      }, options = list(pageLength = 10))
    # }
  })


  ## Cluster Heatmap
  observeEvent(input$dgeHeatButton, {
    if (!is.null(de$markers)) {
      de$dgePlot <-
        getClusterHeatmap(finData$finalData, de$markers, input$clustHeatInput)

      output$dgePlot <- renderPlot({
        de$dgePlot
      })

      updateTabsetPanel(session, "deMainTabSet", selected = "dePlotTab")
    }
  })


  ## DE Plots
  observeEvent(input$dgePlotButton, {
    if (!is.null(de$markers)) {
      de$dgePlot <-
        genePlot(finData$finalData,
                 input$dgePlotType,
                 input$geneNameInput,
                 session)


      output$dgePlot <- renderPlot({
        de$dgePlot
      })

      updateTabsetPanel(session, "deMainTabSet", selected = "dePlotTab")

    } else {
      sendSweetAlert(
        session = session,
        title = "Marker Data Not Found",
        text = "Please run the differential expression pipeline first",
        type = "warning"
      )

    }
  })


  output$downloaddgePlot <- downloadHandler(
    filename = function() {
      paste("DEplot", device = ".png", sep = "")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(
          ...,
          width = width,
          height = height,
          units = "px",
          pointsize = 14
        )
      }
      ggsave(
        file,
        plot = de$dgePlot,
        device = device,
        width = 1280,
        height = 720,
        limitsize = FALSE
      )
    }
  )

  return(de)
}


#' Cluster Heatmap
#'
#' Heatmap generated with Suerat
#'
#' @param s_object Seurat object with clustered data
#' @param markers Differential expression data
#' @param geneNo Number of genes to be displayed
#'
#' @export
#' @return Diffenretial Expression data
getClusterHeatmap <- function(s_object, markers, geneNo) {
  topMarkers <-
    markers %>% group_by(cluster) %>% top_n(n = geneNo, wt = avg_logFC)
  p <- DoHeatmap(s_object, features = topMarkers$gene)
  return(p)
}


#' Gene Plots across clusters
#'
#' Function with exception handling that enables comparison of genes
#' across the different clusters
#'
#' @param finalData Seurat object with cluster data
#' @param plotType The desired plot type
#' @param geneName The name/symbol of the gene of interest
#'
#'
#' @export
#' @return Diffenretial Expression data
genePlot <- function(finalData, plotType, geneName, session) {
  out <- tryCatch(
    {
      if (plotType == 1) {
        VlnPlot(finalData, features = geneName)

      } else if (plotType == 2) {
        FeaturePlot(finalData, features = geneName)

      } else if (plotType == 3) {
        RidgePlot(finalData, features = as.character(geneName))
      }
    },
    error=function(cond) {

      sendSweetAlert(
          session = session,
          title = "Gene not found",
          text = "Please enter an existing gene name/symbol",
          type = "error"
      )

      return(NA)      # Choose a return value in case of error
    }
  )
  return(out)
}
