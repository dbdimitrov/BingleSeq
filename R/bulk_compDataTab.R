#' Bulk Compare Data UI
#'
#' @export
#' @return None
bulk_compDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    sidebarPanel(

      h4("Initialize by Running the DE Pipelines"),

      actionButton(ns("comparisonButton"), "Run DE Pipelines"),

      conditionalPanel(
        condition = "input.comparisonButton > 0",
        ns = ns,

        tags$hr(),

        h4("Generate Venn Diagram"),

        actionButton(ns("vennButton"), "Generate Venn Diagram"),

        conditionalPanel(condition =  "input.vennButton > 0",
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
        ),

        tags$hr(),

        h4("Produce Ranking Consesus"),

        actionButton(ns("rankButton"), "Generate Ranking Table")
      )
    ),

    mainPanel(

      tabsetPanel(
        id = ns("compMainTabSet"),
        tabPanel(
          title = "Venn Diagram",
          value = "compPlotTab",
          htmlOutput(ns("helpCompInfo")),
          plotOutput(ns("comparsionPlot"), width = "800px", height = "500px")
        ),
        tabPanel(title = "Rank Table",
                 value = "compTableTab",
                 DT::dataTableOutput(ns("rankTable")),
                 conditionalPanel(condition = "input.rankButton > 0",
                                  ns = ns,
                                  downloadButton(ns("downloadRank"), "Download Ranking Consensus")
                 )
        )
      )
    )
  )
}


#' Bulk Compare Data Server
#' @param rv filtered counts
#' @param de DE results and meta
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
        the overlap between selected DE methods can be generated <br>
        Each overlap(intersect) can then be downloaded. <br>
        Furthermore, a ranking consesus between the packages can be
        generated and downloaded. <br><br>

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

    rv$allDE <- getAllDE(rv$filteredCounts, de$deTable[[3]], de$batched)

    hide_waiter()
  })


  observeEvent(input$vennButton, {

    rv$xlist <- generateIntersect(rv$allDE)

    output$comparsionPlot <- renderPlot({
      grid.newpage()
      grid.draw(plotAllVenn(rv$xlist))

    })

    updateTabsetPanel(session,
                      "compMainTabSet",
                      selected = "compPlotTab")

  })


  observeEvent(input$rankButton, {

    rv$consensus <- rankConsesus(rv$allDE[[2]][[1]], rv$allDE[[3]][[1]], rv$allDE[[1]][[1]], 2)

    output$rankTable <-
      DT::renderDataTable(rv$consensus)

    updateTabsetPanel(session,
                      "compMainTabSet",
                      selected = "compTableTab")

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

  output$downloadRank <- downloadHandler(
    filename = function() {
      paste("RankingConsesus_",
            ".csv",
            sep = "")
    },
    content = function(file) {
      data <- rv$consensus

      write.csv(data, file, row.names = FALSE)
    }
  )

}



#' Bulk Generate Data required to compare DE Package Results
#'
#' @param readCounts Filtered Counts Table
#' @param meta Metadata table
#' @param batched, whether batch effect was applied to the tables
#' @export
#' @return Returns a list with DE genes according to the different packages
getAllDE <- function(readCounts, meta, batched){
  x1 <- deSequence(readCounts, meta, "Wald", "parametric", batched)
  x2 <- deEdgeR(readCounts, meta, "exact", "TMMwsp", batched)
  x3 <- deLimma(readCounts, meta, "TMMwsp", batched)

  delist <- list(x1, x2, x3)

  return(delist)
}



#' Extract Intersects
#'
#' @param deList list containing the DE Results of the 3 packages
#' @export
#' @return Returns a list with DE genes according to the different packages
generateIntersect <- function(deList) {

  x1_sig <- subset(deList[[1]][[1]], FDR < 0.05)
  x2_sig <- subset(deList[[2]][[1]], FDR < 0.05)
  x3_sig <- subset(deList[[3]][[1]], FDR < 0.05)


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


#' Generates a ranking consesus from the selected methods
#'
#' @param x1_data DE method 1 results
#' @param x2_data DE method 2 results
#' @param x2_data DE method 3 results
#' @param pipeline The pipeline used - scRNA-Seq or Bulk
#' @export
#' @return Returns a dataframe with the ranking consensus
rankConsesus <- function(x1_data, x2_data, x3_data, pipeline){


  if(pipeline == 1){

    #rename p_val_adj to FDR
  }
  x1_data$X <- row.names(x1_data)
  x2_data$X <- row.names(x2_data)
  x3_data$X <- row.names(x3_data)


  # order the data according to most significant
  edgeR_ord <- x1_data[order(x1_data$FDR), ]
  limma_ord <- x2_data[order(x2_data$FDR), ]
  DESeq_ord <- x3_data[order(x3_data$FDR), ]





  x1 <- as.vector(edgeR_ord$FDR) # extract FDR
  names(x1) <- edgeR_ord$X # assign names
  x1_ranked <- rank(x1) # produce ranks
  x1_ranked <- x1_ranked[order(names((x1_ranked)))]

  x2 <- as.vector(limma_ord$FDR)
  names(x2) <- limma_ord$X
  x2_ranked <- rank(x2)
  x2_ranked <- x2_ranked[order(names((x2_ranked)))]

  x3 <- as.vector(DESeq_ord$FDR)
  names(x3) <- DESeq_ord$X
  x3_ranked <- rank(x3)
  x3_ranked <- x3_ranked[order(names((x3_ranked)))]


  # combine the vectors into a
  xdf <- as.data.frame((cbind(x1_ranked,x2_ranked,x3_ranked)))


  # produce consesus
  xdf$consensus <- (rowSums(xdf[,1:3])/3)


  consesus <- as.vector(xdf$consensus) # extract consesus
  names(consesus) <- row.names(xdf) # assign names
  consesus


  rerank <- rank(consesus)
  rerank <- rerank[order(names((rerank)))]

  # p.adj vectors
  edgeR_padj <- as.vector(x1_data$FDR)
  names(edgeR_padj) <- x1_data$X
  edgeR_padj <- edgeR_padj[order(names((edgeR_padj)))]

  limma_padj <- as.vector(x2_data$FDR)
  names(limma_padj) <- x2_data$X
  limma_padj <- limma_padj[order(names((limma_padj)))]

  DESeq_padj <- as.vector(x3_data$FDR)
  names(DESeq_padj) <- x3_data$X
  DESeq_padj <- DESeq_padj[order(names((DESeq_padj)))]


  # rebind the rankings + FDR for each package + consesus
  xdf <- as.data.frame(cbind(x1_ranked, edgeR_padj,
                             x2_ranked, limma_padj,
                             x3_ranked, DESeq_padj,
                             rerank))




  xdf <- xdf[order(xdf$rerank),]

  if(pipeline == 2){
    colnames(xdf) <- c("edgeR Rank", "edgeR adj.p-value",
                       "limma Rank", "limma adj.p-value",
                       "DESeq2 Rank", "DESeq2 adj.p-value",
                       "Ranking Consesus")
  } else {
    colnames(xdf) <- c("T-test Rank", "T-test adj.p-value",
                       "Wilcoxon Rank", "Wilcoxon adj.p-value",
                       "MAST Rank", "MAST adj.p-value",
                       "Ranking Consesus")
  }

  return(xdf)
}
