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

      numericInput(ns("filterValue"), label = "Filter Genes with CPM below", value = 1),

      conditionalPanel(
        condition = "input.selectFilter == 1",
        ns = ns,
        numericInput(ns("sampleNoInput"), label = "in at least # samples", value = 3)

      ),

      fluidRow(column(3, verbatimTextOutput(
        ns("filterValue")
      ))),

      actionButton(ns("filterButton"), label = "Filter Data"),

      tags$hr(),

      conditionalPanel(
        condition = "input.filterButton > 0 && output.panelStatus",
        ns = ns,

        h4("Batch Effect Correction"),

        selectInput(
          ns("selectBatchMethod"),
          label = "Select Method",
          choices = list("Harman" = 1, "ComBat" = 2)
        ),

        actionButton(ns("batchButton"), label = "Correct Batch Effect")
      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(tabsetPanel(
      id = ns("qcTabSet"),

      tabPanel(title = "Filter Data",
          fluidRow(
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
        )
      ),
      tabPanel(title = "Batch Effect Correction",
               value = "batchTab",
               htmlOutput(ns("preBatchHelp")),
               plotOutput(ns("preBatchPCA"), width = "800px", height = "500px"),

               conditionalPanel(
                     condition = "input.filterButton > 0",
                     ns = ns,
                     tags$hr(),

                     htmlOutput(ns("postBatchHelp")),
                     plotOutput(ns("postBatchPCA"), width = "800px", height = "500px")
               )
      )
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
  filt$batchCorrected <- NULL

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


    # batch Effect is on only when more than 1 batch group
    if(length(levels(as.factor(counts$metaTable$batch))) > 1){


      output$preBatchHelp <- renderUI({
        HTML("<h4 style='padding-top: 8px'>Uncorrected Batch Samples PCA</h4>
             <p style='padding-bot: 8px;'><i>
              Samples are coloured by their Batch Groups.</i></p>")
      })

      output$preBatchPCA <- renderPlot({
        plotPCA(filt$filteredCounts, FALSE, counts$metaTable, "batch")
      })

      output$panelStatus <- reactive({
        input$select1=="show"
      })

      outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)
    } else {

      output$panelStatus <- reactive({
        FALSE
      })

      outputOptions(output, "panelStatus", suspendWhenHidden = FALSE)
    }


    # Used to generate DE Tab only when generateSummary is OK
    filt$correctFormat <- TRUE
  })


  observeEvent(input$selectFilter , {
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


  observeEvent (input$batchButton, {
    if (input$selectBatchMethod == 1) {
      filt$batchCorrected <- batchHarman(filt$filteredCounts, counts$metaTable, session)
    } else if (input$selectBatchMethod == 2) {
      filt$batchCorrected <- batchComBat(filt$filteredCounts, counts$metaTable, session)
    }

    if(!is.null(filt$batchCorrected)){

      output$postBatchHelp <- renderUI({
        HTML("<h4 style='padding-top: 8px'>Corrected Batch Samples PCA</h4>
             <p style='padding-bot: 8px;'><i>
              Samples are coloured by their Batch Groups.</i></p>")
      })

      output$postBatchPCA <- renderPlot({
        plotPCA(filt$batchCorrected, TRUE, counts$metaTable, "batch")
      })

      updateTabsetPanel(session,
                        "qcTabSet",
                        selected = "batchTab")
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


#' Harman Batch Effect Correction
#'
#' @param counts, Filtered Count Table
#' @param meta, Metadata Table
#'
#' @export
#' @return Returns batch-effect corrected table
#' Harman Batch Effect Correction
#'
#' @param counts, Filtered Count Table
#' @param meta, Metadata Table
#'
#' @export
#' @return Returns batch-effect corrected table
batchHarman <- function(counts, meta, session){

  out <- tryCatch(
    {
      head(counts)
      countTable  <- counts[, -1]
      rownames(countTable) <- counts[, 1]

      head(countTable)

      harman_results <- harman(countTable, meta$treatment, meta$batch, limit=0.95)

      harman_corr <- reconstructData(harman_results)
      harman_corr[harman_corr<0] <- 0
      out <- harman_corr
    },
    error=function(cond) {

      sendSweetAlert(
        session = session,
        title = "Batch Effect Correction Error",
        text = "Please ensure that the metadata table is in the correct format.",
        type = "error"
      )

      return(NA)
    }
  )

  return(out)
}


#' ComBat Batch Effect Correction
#'
#' @param counts, Filtered Count Table
#' @param meta, Metadata Table
#'
#' @export
#' @return Returns batch-effect corrected table
#' ComBat Batch Effect Correction
#'
#' @param counts, Filtered Count Table
#' @param meta, Metadata Table
#'
#' @export
#' @return Returns batch-effect corrected table
batchComBat <- function(counts, meta, session){

  out <- tryCatch(
    {
      countTable  <- counts[, -1]
      rownames(countTable) <- counts[, 1]

      design <- model.matrix(~as.factor(treatment), data=meta)
      out = ComBat(as.matrix(countTable), batch=meta$batch, mod=design)
    },
    error=function(cond) {

      sendSweetAlert(
        session = session,
        title = "Batch Effect Correction Error",
        text = "Please ensure that the metadata table is in the correct format.",
        type = "error"
      )

      return(NA)
    }
  )

  return(out)
}



#' Plot PCA
#'
#' Uses factoextra package to plot a PCA plot
#'
#' @param data Differential Expression results (deTable)
#' @param rowNames, Boolean checking whether genes are row names or not
#' @param meta, metadata table
#' @param col, colour by batch or treatment
#' @export
#' @return Returns a PCA plot
plotPCA <- function(data, rowNames, meta, col) {

  if(!rowNames){
    countTable  <- data[, -1]
    rownames(countTable) <- data[, 1]
  } else {
    countTable <- data
  }

  if(is.data.frame(data)){
    countTable <- as.matrix(sapply(countTable, as.numeric))
  }

  xx <- prcomp(t(countTable))

  if(col == "batch"){
    group <- meta$batch
    leg <- "Batch"
  } else if(col == "treatment"){
    group <- meta$treatment
    leg <- "Treatment"
  }

  return(fviz_pca_ind(xx,
                      repel = FALSE,
                      habillage=group,
                      title = "",
                      legend.title = leg,
                      axes.linetype=NA) +
           theme_classic(base_size = 16))
}


#' Plot Scree
#'
#' Uses factoextra package to plot a Scree plot
#'
#' @param data Differential Expression results (deTable)
#' @export
#' @return Returns a Scree plot
plotScree <- function(data) {

  x <- as.matrix(sapply(data, as.numeric))
  xx <- prcomp(t(x))

  return(fviz_eig(xx,
                  title=""
                  ) + theme_classic(base_size = 16))
}
