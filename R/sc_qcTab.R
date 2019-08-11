#' Single Cell Quality Control UI
#'
#' @export
#' @return None
sc_qcUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(
      h4("Prefilter Cells and Counts"),

      numericInput(
        ns("minCellsObject"),
        label = "Minimum number of cells per gene",
        min = 1,
        value = 3
      ),
      numericInput(
        ns("minGenesObject"),
        label = "Minimum number of genes per cell",
        min = 1,
        value = 200
      ),


      actionButton(ns("preqcButton"), label = "Initialize Project"),

      conditionalPanel(
        condition =  "input.preqcButton > 0",
        ns = ns,

        tags$hr(),

        h4("Filter Cell outliers"),

        numericInput(
          ns("minFeatureInput"),
          label = "Minimum Feature No",
          min = 1,
          value = 100
        ),



        numericInput(
          ns("maxFeatureInput"),
          label = "Maximum Feature No",
          min = 1,
          value = 5000
        ),

        actionButton(ns("postqcButton"), label = "Filter Cells")

      )

    ),

    # Main panel for displaying outputs ----
    mainPanel(tabsetPanel(
      id = ns("qcTabSet"),

      #Add Text with No of cells and features as lodaded

      tabPanel(title = "Object Preview",

               htmlOutput(ns("helpQCInfo")),

               DT::dataTableOutput(ns("dataTable"))
               ),


      tabPanel(
        title = "Outlier Violin Plot",
        value = "tab2_val",
        verbatimTextOutput(ns("preFilterText"), placeholder = T),
        plotOutput(ns("preqcPlot"), width = "100%", height = "500px"),

        tags$hr(),

        verbatimTextOutput(ns("postFilterText"), placeholder = T),
        plotOutput(ns("postqcPlot"), width = "100%", height = "500px")
      )
    ))
  )
}

#' Single Cell Quality Control UI
#'
#' @param countsT countTable loaded by loadTab
#' @export
#' @return A reactive value contaning the filtered data
sc_qc <- function(input, output, session, countsT) {
  filt <- reactiveValues()


  output$helpQCInfo <- renderUI({
    if(is.null(filt$data)){
      HTML("<div style='border:2px solid blue; padding-top: 8px; padding-bot: 8px; font-size: 14px;
      border-radius: 10px;'>
      <p style='text-align: center'><b>This tab enables Quality control. </b> </p> <br>
      First, filter genes detected below a certain number of cells and cells with less than a certain number of expressed genes. <br>
      <i> This is done as the project (object) is initialized and subsequently displayed as a table. </i> <br>
      Then, swap to the 'Outlier Violin plot' tab to visualize and exclude Cell Outliers. </div>")
    } else{
      HTML("")
    }

  })


  ### Show PreQC INFO -------
  observeEvent(input$preqcButton, {
    filt$data <- CreateSeuratObject(
      counts = countsT$countTable,
      min.cells = input$minCellsObject,
      min.features = input$minGenesObject,
      project = ""
    )


    output$dataTable <- DT::renderDataTable(DT::datatable(as.data.frame(filt$data@meta.data)))

    output$preFilterText <- renderText({
      sprintf("Gene(Feature) No: %s;   Cell No: %s.",
              nrow(filt$data),
              ncol(filt$data))
    })

    output$preqcPlot <- renderPlot({
      VlnPlot(
        filt$data,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2
      )

    })
  })



  ### SHOW postQC INFO ------
  observeEvent(input$postqcButton, {
    filt$filteredData <- seuratQC(filt$data, input$minFeatureInput,
                                  input$maxFeatureInput, session)


    output$postFilterText <- renderText({
      sprintf(
        "Gene(Feature) No: %s;   Cell No: %s.",
        nrow(filt$filteredData),
        ncol(filt$filteredData)
      )
    })

    output$postqcPlot <- renderPlot({
      VlnPlot(
        filt$filteredData,
        features = c("nFeature_RNA", "nCount_RNA"),
        ncol = 2
      )
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
#' @param session Current R session supplied by the server
#'
#' @export
#' @return A Seurat_object with the filtered data
seuratQC <- function(seurat_object, minF, maxF, session) {


  seurat_object <- tryCatch(
    {
      keep <- seurat_object$nFeature_RNA > minF &
        seurat_object$nFeature_RNA < maxF


      seurat_object <- seurat_object[, keep]
    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data format error",
        text = "Ensure that a correctly formatted data was supplied",
        type = "error"
      )
      return()
    }
  )

  return(seurat_object)
}
