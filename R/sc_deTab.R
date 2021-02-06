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

              h4("Generate DE Data"),

              selectInput(
                ns("dgeTestCombo"),
                label = "Select Test Type",
                choices = list(
                  "MAST" = "MAST",
                  "Wilcoxon Rank Sum test" = "wilcox",
                  "Student's T-test" = "t",
                  "DESeq2" = "DESeq2",
                  "Logistic Regression" = "LR"
                )
              ),

              numericInput(
                ns("logFCdgeInput"),
                label = "Fold-Change Threshold >",
                min = 0,
                max = 10,
                value = 2
              ),

              numericInput(
                ns("adjPdgeInput"),
                label = "Adjusted P-value Threshold <",
                min = 0,
                max = 1,
                value = 0.05
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
                     htmlOutput(ns("helpDEInfo")),
                     DT::dataTableOutput(ns("dgeTable")),
                     conditionalPanel(condition = "input.dgeButton > 0",
                                      ns = ns,

                                      downloadButton(ns("deDownload"), "Download Table")

                     )

                     ),
            tabPanel(
              title = "Plot",
              value = "dePlotTab",

              plotOutput(ns("dgePlot"), width = "1280px", height = "720px"),
              downloadButton(ns("downloaddgePlot"), "Download Plot")
            )
          )))
}

#' Single Cell Differential Expession Tab Server
#'
#' @param finData Reactive value containing a seurat object with clustered data
#'
#' @export
#' @return  Returns Diffenretial Expression data
sc_de <- function(input, output, session, finData) {
  de <- reactiveValues()

  output$helpDEInfo <- renderUI({
    if(input$dgeButton == 0){
      HTML(
        "<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px;'>

        <p style='text-align: center'>
        <b>This tab enables DE analysis of clustering results.</b> </p> <br>
        To indentify Marker Genes for each cluster,
        proceed first by selecting the preferred DE method. <br>
        Then specify pre-filter options according to: <br>
        Fold-change, adj. P-value threshold,
        and genes expressed in a minimum fraction of cells. <br> <br>
        <i>Note: MAST was shown to be among
        the best scRNA-Seq DE methods (Soneson & Robinson, 2018),
        as such it is likely the best option here. </i>
        Also, please be patient as DE analysis is run on all clusters
        and as such it may take some time.
        MAST is particularly time consuming.  </div>"
        
      )
    } else {
      HTML("")
    }
  })

  ## Generate DE Data
  observeEvent(input$dgeButton, {
    waiter_show(html=tagList(spin_folding_cube(), h2("Loading...Stay Patient :)")))

    if(input$dgeTestCombo == "DESeq2"){
      finData$finalData[["RNA"]]@counts <- as.matrix(finData$finalData[["RNA"]]@counts) + 1
    }

    de$markers <- FindAllMarkers(
      finData$finalData,
      test.use = input$dgeTestCombo,
      min.pct = input$pctdgeInput,
      logfc.threshold = log(input$logFCdgeInput)
    )

    # Filter by adjusted P-value
    filter <- de$markers$p_val_adj < input$adjPdgeInput
    de$markers <- de$markers[filter,]


    waiter_hide()

    write.csv(de$markers, file=paste0(tempdir(),
                                      "/AllMarkerGenes_",
                                      input$dgeTestCombo,
                                      ".csv"), row.names = FALSE)

    output$dgeTable <-
      DT::renderDataTable(if (input$dgeClusterCheck) {
        de$markers[,1:(ncol(de$markers)-1)] %>% 
          rownames_to_column("gene_id") %>% 
          datatable(rownames = FALSE)
      } else{
        de$markers[de$markers$cluster == input$dgeClustInput,
                   1:(ncol(de$markers)-1)] %>%
          rownames_to_column("gene_id") %>% 
          datatable(rownames = FALSE)
      }, options = list(pageLength = 10))
  })


  ## Cluster Heatmap
  observeEvent(input$dgeHeatButton, {
    if (!is.null(de$markers)) {
      waiter_show(html=tagList(spin_folding_cube(), h2("Loading ...")))

      de$dgePlot <-
        getClusterHeatmap(finData$finalData, de$markers, input$clustHeatInput)

      hm.palette <-
        colorRampPalette(c("red", "white", "blue")) # Set the colour range

      de$dgePlot <- de$dgePlot +
        scale_fill_gradientn(colours = hm.palette(100))
      

      output$dgePlot <- renderPlot({
        de$dgePlot
      })

      waiter_hide()

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
        de$dgePlot  +
          theme(axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),  
                axis.title.x = element_text(size = 16),
                axis.title.y = element_text(size = 16))
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

  output$deDownload <- downloadHandler(
    filename = function() {
      paste(format(Sys.time(), "%y-%m-%d_%H-%M"), "_deResults" , ".csv", sep = "")
    },
    content = function(file) {
      data <- de$markers

      write.csv(data, file)
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
#' @return Returns DE Heatmap
getClusterHeatmap <- function(s_object, markers, geneNo) {
  topMarkers <-
    markers %>% 
    group_by(cluster) %>%
    top_n(n = geneNo, wt = avg_logFC)
  p <- DoHeatmap(s_object, features = topMarkers$gene)  +
    theme_classic(base_size = 14)
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
#' @return Returns DE Gene Plots
genePlot <- function(finalData, plotType, geneName, session) {
  out <- tryCatch(
    {
      if (plotType == 1) {
        VlnPlot(finalData, features = geneName) +
          theme_classic(base_size = 20)

      } else if (plotType == 2) {
        FeaturePlot(finalData, features = geneName, pt.size = 1.5) +
          theme_classic(base_size = 20)

      } else if (plotType == 3) {
        RidgePlot(finalData, features = as.character(geneName))  +
          theme_classic(base_size = 20)
      }
    },
    error=function(cond) {

      sendSweetAlert(
          session = session,
          title = "Gene not found",
          text = "Please enter an existing gene name/symbol",
          type = "error"
      )

      return(NA)
    }
  )
  return(out)
}
