#' Bulk Compare Data UI
#'
#' @export
#' @return None
bulk_compDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarPanel(
      h4("Generate Venn Diagram"),

      actionButton(ns("comparisonButton"), "Run DE Pipelines"),

      conditionalPanel(
        condition = "input.comparisonButton > 0",
        ns = ns,

        selectInput(
          ns("selectIntersect"),
          label = ("Select intersect of interest"),
          choices = list(
            "DESeq2 (only)" = "a1",
            "DESeq2 & edgeR" = "a2",
            "edgeR (only)" = "a3",
            "DESeq2 & limma" = "a4",
            "All Three Packages" = "a5",
            "edgeR & limma" = "a6",
            "limma (only)" = "a7"
          )
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


#' Bulk Compare Data Server
#' @param rv filtered counts
#' @param de DE results and replicates/conditions info
#' @export
#' @return None
bulk_compData <- function(input, output, session, rv, de) {
  output$helpCompInfo <- renderUI({
    if (input$comparisonButton == 0) {
      HTML(
        "<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px;'>

        <p style ='text-align: center'> <b> This tab supplies
        users with an option to assess the agreement
        between the different DE analysis packages.</b> </p> <br>

        Once the pipeline is finished a Venn Diagram with
        the overlap between selected DE methods is returned. <br>
        Each overlap(intersect) can then be downloaded <br> <br>

        <i>Note that Differentially Expressed Genes are considered
        significant if FDR adjusted p-value < 0.05. <br>
           Moreover, the procedure runs each DE analysis pipeline,
        as such it is rather time-consuming.</i> </div>"
      )
    } else {
      HTML("")
    }
  })

  observeEvent(input$comparisonButton, {
    show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))

    rv$xlist <-
      getAllDE(rv$filteredCounts, de$conditionNo, de$replicateNo)

    output$comparsionPlot <- renderPlot({
      grid.newpage()
      grid.draw(plotAllVenn(rv$xlist))

    })

    hide_waiter()
  })


  observeEvent(input$selectIntersect, {
    if (!is.null(rv$xlist)) {
      rv$intersect <- getIntersect(rv$xlist, input$selectIntersect)
    }
  })

  output$downloadIntersect <- downloadHandler(
    filename = function() {
      paste("PackageComparison_intersect_",
            input$selectIntersect,
            ".csv",
            sep = "")
    },
    content = function(file) {
      data <- rv$intersect

      write.csv(data, file, row.names = FALSE)
    }
  )
}



#' Bulk Generate Data required to compare DE Package Results
#'
#' @param readCounts Filtered Counts Table
#' @param conditionNo Number of Conditions
#' @param replicateNo Number of Replicates
#' @export
#' @return Returns a list with DE genes according to the different packages
getAllDE <- function(readCounts, conditionNo, replicateNo) {
  x1 <- deSequence(readCounts, conditionNo, replicateNo)
  x2 <- deEdgeR(readCounts, conditionNo, replicateNo)
  x3 <- deLimma(readCounts, conditionNo, replicateNo)

  x1_sig <- subset(x1, FDR < 0.05)
  x2_sig <- subset(x2, FDR < 0.05)
  x3_sig <- subset(x3, FDR < 0.05)


  DESeq_genes <- as.vector(rownames(x1_sig))
  EdgeR_genes <- as.vector(rownames(x2_sig))
  Limma_genes <- as.vector(rownames(x3_sig))

  list <- list(DESeq_genes, EdgeR_genes, Limma_genes)

  names(list) <- c("DESeq2", "edgeR", "limma")

  return(list)
}

#' Generate a Venn Diagram with DE genes
#'
#' Visual representation of Package Agreement/Disagreement
#'
#' @param xlist A list with DE genes according to the different packages
#' @export
#' @return Returns a Venn Diagram
plotAllVenn <- function(xlist) {
  venn.plot <- venn.diagram(
    xlist,
    filename = NULL,
    fill = c("red", "yellow", "blue"),
    lty = "blank",
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.6,
    cat.fontface = "bold",
    cat.default.pos = "outer"
  )


  return(venn.plot)
}

#' Returns the Genes of a given intersect within the Venn Diagram
#'
#' @param xlist A list with DE genes according to the different packages
#' @param intersectID The corresponding number of the intersect of interest
#' @export
#' @return Returns the names of the genes within that intersect
getIntersect <- function(xlist, intersectID) {
  x <- calculate.overlap(xlist)

  inter <- x[[intersectID]]

  return(inter)
}
