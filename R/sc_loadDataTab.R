#' Single Cell Load Data UI
#'
#' @export
#' @return None
sc_loadDataUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Sidebar layout with input and output definitions ----
    sidebarPanel(
      
      h4("Upload Data"),
      
      radioButtons(
        ns("loadData"),
        label = "Data Type",
        c("Count Data Table" = 1,
          "10x Genomics Data" = 2,
          "Test Single cell Data" = 3)
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
        textInput(ns("location10xInput"),
                  "10x folder location:",
                  "For example: ~/10x_test"),
        actionButton(ns("directoryButton"), label = "Load 10x Data")
      ),
      
      conditionalPanel(
        condition = "input.loadData == 3",
        ns = ns,
        tags$hr(),
        actionButton(ns("scTestDataButton"), label = "Load Test Data")
        
      )
      
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(
      htmlOutput(ns("helpLoadInfo")),
      verbatimTextOutput(ns("loadDataText"), placeholder = F)
    )
  )
}


#' Single Cell Load Data Server
#'
#' Enables the upload of two types of data: 10x Genomics and Count data
#'
#' @export
#' @return Returns a count table to be used in the sc analysis
sc_loadData <- function(input, output, session) {
  counts <- reactiveValues()
  
  
  
  output$helpLoadInfo <- renderUI({
    if(as.numeric(input$loadData)==1 && is.null(counts$countTable)){
      HTML("<div style='border:2px solid blue; font-size: 14px;
      border-radius: 10px;text-align: center'>
      <p style='padding-top: 8px'> Select a count table
      that contains the read counts in .csv/.txt file format. </p>
      <p style ='font-style: italic; padding-bottom: 8px;'> 
      Note: The first row(header) and first column should contain cell names 
      and gene names/IDs, respsectively </p> </div>")
    }else if(as.numeric(input$loadData)==2 && is.null(counts$countTable)){
      HTML("<div style='border:2px solid blue; font-size: 14px;
      border-radius: 10px;text-align: center'>
           <p style='padding-top: 8px'; padding-bottom: 8px;>
           Select a 10x Genomics output directory containing
           matrix.mtx, barcodes.tsv, and genes.tsv files </p> </div>")
    }else if(as.numeric(input$loadData)==3 && is.null(counts$countTable)){
      HTML("<div style='border:2px solid blue;
      font-size: 14px;
      border-radius: 10px;text-align: center'>
           <p style='padding-top: 8px'; padding-bottom: 8px;>
           Load an example public Cellranger 10x Genomics dataset with 3k PBMCs
           </p> </div>")
    } else{
      HTML("")
    }
  })
  
  # Load 10x -----
  observeEvent(input$directoryButton, {
      print(input$location10xInput)
      waiter_show(tagList(spin_folding_cube(), h2("Loading ...")))
      counts$countTable <- load10xData(input$location10xInput, session)
      hide_waiter()
  })
  
  # Load CountTable -----
  observeEvent(input$file1, {
    waiter_show(tagList(spin_folding_cube(), h2("Loading ...")))
    counts$countTable <- read.csv(input$file1$datapath,
                                  sep = input$sep,
                                  row.names = 1)
    
    hide_waiter()
    
  })
  
  
  # Load Test Data -----
  observeEvent(input$scTestDataButton, {
    sc_example_data <-paste0(system.file("extdata",
                                         "hg19",
                                         package = "BingleSeq"), "/") 
    req(nchar(sc_example_data > 0))
    counts$countTable <- load10xData(sc_example_data, session)
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
#' @return Returns 10x Genomics Data
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
