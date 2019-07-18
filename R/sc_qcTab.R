#' Single Cell Quality Control UI
#'
#' @export
#' @return None
sc_qcUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(

      h5("Prefilter Features:"),

      numericInput(ns("minCellsObject"), label = "Mininum cells expressed",
                   min = 1, value = 3),
      numericInput(ns("minGenesObject"), label = "Minimum summed counts",
                   min = 1, value = 200),


      actionButton(ns("preqcButton"), label = "Initialize Project"),

      tags$hr(),

      shinyjs::hidden(


            numericInput(ns("minFeatureInput"), label = "Minimum Feature No",
                         min = 1, value = 100),



            numericInput(ns("maxFeatureInput"), label = "Maximum Feature No",
                         min = 1, value = 5000),


            actionButton(ns("postqcButton"), label = "Filter Data")

      )

    ),

    # Main panel for displaying outputs ----
    mainPanel(

      tabsetPanel(id = ns("qcTabSet"),

                  #Add Text with No of cells and features as lodaded

                  tabPanel(title = "Object Preview",
                           DT::dataTableOutput(ns("dataTable"))
                  ),


                  # extendShinyjs(script = "javascript.js"),

                  tabPanel(title = "QC Violin Plot", value = "tab2_val",
                           verbatimTextOutput(ns("preFilterText"), placeholder = T),
                           plotOutput(ns("preqcPlot")),

                           tags$hr(),

                           verbatimTextOutput(ns("postFilterText"), placeholder = T),
                           plotOutput(ns("postqcPlot"))
                  )
      )


    )
  )
}

#' Single Cell Quality Control UI
#'
#' @param countsT countTable loaded by loadTab
#' @export
#' @return A reactive value contaning the filtered data
sc_qc <- function(input, output, session, countsT) {

  filt <- reactiveValues()

  ### Show PreQC INFO -------
  observeEvent(input$preqcButton, {

    filt$data <- CreateSeuratObject(counts = countsT$countTable,
                                    min.cells = input$minCellsObject,
                                    min.features = input$minGenesObject,
                                    project = "userProject")


    output$dataTable <- DT::renderDataTable(
      DT::datatable(as.data.frame(filt$data@meta.data))

    )

    output$preFilterText <- renderText({

      sprintf("Gene(Feature) No: %s;   Cell No: %s.", nrow(filt$data), ncol(filt$data))
    })

    output$preqcPlot <- renderPlot({

      VlnPlot(filt$data, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

    })

    shinyjs::showElement("minFeatureInput")
    shinyjs::showElement("maxFeatureInput")
    shinyjs::showElement("postqcButton")

  })



  ### SHOW postQC INFO ------
  observeEvent(input$postqcButton, {

    filt$filteredData <- seuratQC(filt$data, input$minFeatureInput,
                                  input$maxFeatureInput)


    output$postFilterText <- renderText({
      sprintf("Gene(Feature) No: %s;   Cell No: %s.", nrow(filt$filteredData), ncol(filt$filteredData))
    })

    output$postqcPlot <- renderPlot({

      VlnPlot(filt$filteredData, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
    })
  })

  return(filt)
}


#' Single Cell Quality Control
#'
#' Performs QC on Seurat object, which serves as a backbone for the sc pipeline
#'
#' @param seurat_object The initialized seurat object with the count data
#' @param minF Minimum feature counts per cell
#' @param maxF Maximum feature counts per cell
#' @export
#' @return A Seurat_object with the filtered data
seuratQC <- function(seurat_object, minF, maxF){

  keep <- seurat_object$nFeature_RNA > minF &
    seurat_object$nFeature_RNA < maxF


  seurat_object <- seurat_object[,keep]


  return(seurat_object)
}
