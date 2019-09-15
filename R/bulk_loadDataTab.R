#' Bulk load Data UI
#'
#' @export
#' @return None
bulk_loadDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar layout with input and output definitions ----
    sidebarPanel(

      h4("Load Counts Data"),

      # Input: Select a file ----
      fileInput(
        ns("file1"),
        "Choose Counts File",
        multiple = FALSE,
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      ),

      checkboxInput(ns("header"), tags$b("Header"), FALSE),

      radioButtons(
        ns("sep"),
        "Separator",
        choices = c(
          Comma = ",",
          Space = " ",
          Tab = "\t",
          Semicolumn = ";"
        ),

        selected = ","
      ),

      tags$hr(),

      ######################
      h4("Load Meta Data"),

      # Input: Select a metadata file ----
      fileInput(
        ns("metaFile1"),
        "Choose Counts File",
        multiple = FALSE,
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      ),

      radioButtons(
        ns("metaSep"),
        "Separator",
        choices = c(
          Comma = ",",
          Space = " ",
          Tab = "\t",
          Semicolumn = ";"
        ),

        selected = ","
      )


    ),

    # Main panel for displaying outputs ----
    mainPanel(
              htmlOutput(ns("helpLoadInfo")),
              DT::dataTableOutput(ns("dto")),
              tags$br(),
              htmlOutput(ns("metaInfo")),
              DT::dataTableOutput(ns("metaTable")))
  )
}


#' Bulk load Data SERVER
#'
#'
#' @export
#' @return counts the loaded Count Table
bulk_loadData <- function(input, output, session) {
  counts <- reactiveValues()


  output$helpLoadInfo <- renderUI({
    if(is.null(counts$countTable)){
      HTML("<div style='border:2px solid blue; font-size: 14px; border-radius: 10px;'>
      <p style ='font-size: 15px; text-align: center'; padding-top: 8px;>
      <b>Upload Data tab. </b> </p>

      <p>Please select a count table that
      contains the read counts in .csv/.txt file format. </p>
      <p style ='font-style: italic; padding-bottom: 8px;'>
      Note: The first row (header) may contain contain condition(sample) names,
      and first column should contain gene names/IDs. </p> </div>")
    } else{
      HTML("<h4 style='padding-top: 8px'>Preview Count Table</h4>
           <p style='padding-bot: 8px;'><i>Please ensure
           that the table was loaded appropraitely.</i></p>")
    }
  })



  # Load Counts File ----
  observeEvent(input$file1, {
    counts$countTable <- read.csv(input$file1$datapath,
                                  header = input$header,
                                  sep = input$sep)

    output$dto <- DT::renderDataTable(DT::datatable(counts$countTable, options = list(pageLength = 10)))

  })

  observeEvent(input$metaFile1, {

    if(!is.null(counts$countTable)){
      output$metaInfo <- renderUI({
          HTML("<h4 style='padding-top: 8px'>Preview Meta Data</h4>
             <p style='padding-bot: 8px;'><i>Please ensure
             that the table is formatted appropraitely.</i></p>")
        })
    }

    counts$metaTable <- read.csv(input$metaFile1$datapath,
                                  header = TRUE,
                                  sep = input$metaSep)

    output$metaTable <- DT::renderDataTable(DT::datatable(counts$metaTable, options = list(pageLength = 10)))

  })

  return(counts)

}


#' Generate Summary
#'
#' Generates a summary of the columns(samples) in the Count Table
#'
#' @export
#' @return df a dataframe containing the summary
generateSummary <- function(counts, session) {
  countTable  <- counts[, -1]
  rownames(countTable) <- counts[, 1]

  x1 <- vector()

  x2 <- vector()

  out <- tryCatch(
    {
      for (i in 1:ncol(countTable)) {
        x1[i] <-  colnames(countTable)[i]

        x2[i] <- colSums(countTable[i])
      }

      x1[ncol(countTable) + 1] <- "Total Counts"
      x2[ncol(countTable) + 1] <- sum(x2)

      x1[ncol(countTable) + 2] <- "Sample Median"
      x2[ncol(countTable) + 2] <- median(x2)

      x1[ncol(countTable) + 3] <- "Gene#"
      x2[ncol(countTable) + 3] <- nrow(countTable)


      df <- data.frame(x1, x2)

      format.data.frame(df, big.mark = ",")

      colnames(df) <- c("", "Counts")

      out <- df

    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data Format Error",
        text = "Ensure that correctly formatted data with
        appropriately chosen number of conditions and replicates were supplied",
        type = "error"
      )
      return()
    }
  )
  return(out)
}
