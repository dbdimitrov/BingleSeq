#' Bulk Compare Data UI
#'
#' @export
#' @return None
bulk_compDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarPanel(

      actionButton(ns("comparisonButton"), "Generate Venn Diagram"),



      conditionalPanel(condition = "input.comparisonButton > 0",
                       ns = ns,

                       selectInput(ns("selectIntersect"),
                                   label =("Select intersect of interest"),
                                   choices = list("ALDEx2" = "a1",
                                                  "ALDEx2 & Limma" = "a2",
                                                  "Limma" = "a3",
                                                  "ALDEx2 & Deseq2" = "a4",
                                                  "Limma & ALDEx2 & DESeq2" = "a5",
                                                  "All Four Statistics" = "a6",
                                                  "ALDEx2 & Limma & EdgeR" = "a7",
                                                  "Limma & EdgeR" = "a8",
                                                  "DESeq2" = "a9",
                                                  "DESeq2 & Limma" = "a10",
                                                  "DESeq2 & Limma & EdgeR" = "a11",
                                                  "DESeq2 & ALDEx2 & EdgeR" = "a12",
                                                  "ALDEx2 & EdgeR" = "a13",
                                                  "EdgeR" = "a14",
                                                  "DESeq2 & EdgeR" = "a15")),

                       downloadButton(ns("downloadIntersect"), "Download Intersect")
      )
    ),

    mainPanel(
      plotOutput(ns("comparsionPlot"))
    )
  )
}


#' Bulk Compare Data Server
#' @param rv filtered counts
#' @param de DE results and replicates/conditions info
#' @export
#' @return None
bulk_compData <- function(input, output, session, rv, de) {

  observeEvent(input$comparisonButton, {
    output$comparsionPlot <- renderPlot({

      rv$xlist <- getAllDE(rv$filteredCounts, de$conditionNo, de$replicateNo)

      grid.newpage()
      grid.draw(plotAllVenn(rv$xlist))

    })
  })


  observeEvent(input$selectIntersect, {

    if(!is.null(rv$xlist)){
      rv$intersect <- getIntersect(rv$xlist, input$selectIntersect)
    }


  })

  output$downloadIntersect <- downloadHandler(
    filename = function() {
      paste("PackageComparison_intersect_",input$selectIntersect, ".csv", sep="")
    },
    content = function(file) {

      data <- rv$intersect;
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
#' @return A list with DE genes according to the different packages
getAllDE <- function(readCounts, conditionNo, replicateNo){

  x1 <- deSequence(readCounts, conditionNo, replicateNo)
  x2 <- deEdgeR(readCounts, conditionNo, replicateNo)
  x3 <- deALDE(readCounts, conditionNo, replicateNo)
  x4 <- deLimma(readCounts, conditionNo, replicateNo)

  # with foldchange
  # if(conditionNo==3){
  #   conditionNo = c(2,3,4)
  # }
  #
  #
  # x1_sig <- subset(x1, FDR < 0.05 & abs(x1[conditionNo]) > 1)
  # x2_sig <- subset(x2, FDR < 0.05 & abs(x2[conditionNo]) > 1)
  # x3_sig <- subset(x3, FDR < 0.05 & abs(x3[conditionNo]) > 1)
  # x4_sig <- subset(x4, FDR < 0.05 & abs(x4[conditionNo]) > 1)



  #without

  x1_sig <- subset(x1, FDR < 0.05)
  x2_sig <- subset(x2, FDR < 0.05)
  x3_sig <- subset(x3, FDR < 0.05)
  x4_sig <- subset(x4, FDR < 0.05)


  DESeq_genes <- as.vector(rownames(x1_sig))
  EdgeR_genes <- as.vector(rownames(x2_sig))
  ALDE_genes <- as.vector(rownames(x3_sig))
  Limma_genes <- as.vector(rownames(x4_sig))

  list <- list(DESeq_genes, EdgeR_genes, ALDE_genes, Limma_genes)

  print(list)

  names(list) <- c("DESeq2", "EdgeR", "ALDEx2", "Limma")

  return(list)
}

#' Generate a Venn Diagram with DE genes
#'
#' Visual representation of Package Agreement/Disagreement
#'
#' @param xlist A list with DE genes according to the different packages
#' @export
#' @return A Venn Diagram
plotAllVenn <- function(xlist){

  venn.plot <-venn.diagram(xlist, filename = NULL,
                           fill = c("red", "yellow", "blue", "green"),
                           lty = "blank",
                           fontface = "bold",
                           fontfamily = "sans",
                           cat.cex = 1.5,
                           cat.fontface = "bold",
                           cat.default.pos = "outer")


  return(venn.plot)
}

#' Returns the Genes of a given intersect within the Venn Diagram
#'
#' @param xlist A list with DE genes according to the different packages
#' @param intersectID The corresponding number of the intersect of interest
#' @export
#' @return The names of the genes within that intersect
getIntersect <- function(xlist, intersectID){
  x <- calculate.overlap(xlist)

  inter <- x[[intersectID]]

  return(inter)
}
