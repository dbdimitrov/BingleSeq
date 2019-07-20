#' Bulk load Data UI
#'
#' @export
#' @return None
bulk_loadDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar layout with input and output definitions ----
    sidebarPanel(
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

      tags$hr(),

      # Input: Checkbox if file has header ----
      checkboxInput(ns("header"), "Header", FALSE),

      # Input: Select separator ----
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
      )

    ),

    # Main panel for displaying outputs ----
    mainPanel(tags$h4("Preview Data"),

              DT::dataTableOutput(ns("dto")))
  )
}


#' Bulk load Data SERVER
#'
#'
#' @export
#' @return counts the loaded Count Table
bulk_loadData <- function(input, output, session) {
  counts <- reactiveValues()

  # Load File ----
  observeEvent(input$file1, {
    counts$countTable <- read.csv(input$file1$datapath,
                                  header = input$header,
                                  sep = input$sep)


    output$dto <- DT::renderDataTable(DT::datatable(counts$countTable, options = list(pageLength = 10)))

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

      x1[ncol(countTable) + 1] <- "Total"
      x2[ncol(countTable) + 1] <- sum(x2)

      x1[ncol(countTable) + 2] <- "Median"
      x2[ncol(countTable) + 2] <- median(x2)

      x1[ncol(countTable) + 3] <- "Genes#"
      x2[ncol(countTable) + 3] <- nrow(countTable)


      df <- data.frame(x1, x2)

      format.data.frame(df, big.mark = ",")

      colnames(df) <- c("Sample", "Counts")

      out <- df

      # write.csv(df, file="output/Data_Summary.csv", row.names = FALSE)

    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data Format Error",
        text = "Ensure that correctly formatted data with appropriately chosen number of conditions and replicates were supplied",
        type = "error"
      )
      return()
    }
  )
  return(out)
}
