#' Bulk Plot Data UI
#'
#' @export
#' @return None
bulk_plotDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(

      h4("Visualize DE Results"),

      selectInput(
        ns("selectPlotCombo"),
        label = "Select Plot Type",
        choices = list(
          "PCA plot" = "pca",
          "Scree plot" = "scree",
          "Barchart" = "bar",
          "Volcano Plot" = "volcano",
          "MA Plot" = "MA",
          "Heatmap" = "heat"
        )
      ),

      #VP or MA
      conditionalPanel(
        condition = "input.selectPlotCombo == 'volcano' ||
                               input.selectPlotCombo == 'MA' ||
                               input.selectPlotCombo == 'bar' ||
                               input.selectPlotCombo == 'heat'" ,
        ns = ns,

        numericInput(
          ns("plotPvalue"),
          label = ("Adjusted P-value threshold <"),
          value = 0.05,
          min = 0.0001,
          max = 0.5
        ),
        fluidRow(column(3, verbatimTextOutput(ns(
          "plotPvalue"
        )))),

        numericInput(
          ns("plotFC"),
          label = ("Fold-Change threshold >"),
          value = 2,
          min = 0,
          max = 20
        ),

        verbatimTextOutput(ns("plotFC"))

      ),

      conditionalPanel( # Heatmap
        condition = "input.selectPlotCombo == 'heat'",
        ns = ns,

        numericInput(
          ns("topGeneNo"),
          label = ("Number of genes to display"),
          value = 30,
          min = 10
        ),
        fluidRow(column(3, verbatimTextOutput(ns(
          "topGeneNo"
        )))),

        selectInput(
          ns("selectTypeHM"),
          label = ("Genes of interest:"),
          choices = list(
            "All Significant" = 1,
            "Upregulated" = 2,
            "Downregulated" = 3,
            "Non-significant" = 4
          ),
          selected = 1
        ),


        checkboxInput(ns("topNames"), label = "Include Gene Names",
                      value = TRUE)
      ),


      actionButton(ns("figButton"), "Generate Plot"),

      tags$hr(),

      downloadButton(ns("downloadTopHeat"), "Download Plot"),

      conditionalPanel(
        condition = "input.selectPlotCombo == 'heat'",
        ns = ns,

        downloadButton(ns("downloadTopGenes"), "Download CSV")
      )

    ),


    # Main panel for displaying outputs ----
    mainPanel(
      htmlOutput(ns("helpPlotInfo")),

      plotOutput(ns("mainPlot"), width = "800px", height = "500px"),

      conditionalPanel(
        condition = "input.selectPlotCombo == 'bar'",
        ns = ns
        # ,
        # DT::dataTableOutput(ns("barTable"))
      )
    )
  )
}


#' Bulk Plot Data Tab Server
#'
#' @param rv Reactive value containing the DE results (deTable)
#'
#' @export
#' @return None
bulk_plotData <- function(input, output, session, rv) {

  output$helpPlotInfo <- renderUI({
    if(is.null(rv$plot)){
      HTML("<div style='border:2px solid blue; font-size: 14px; border-radius: 10px;'>

      <p style='text-align: center'><b> This tab enables the
      visualization of DE results. </b> </p> <br>

      <b>PCA plot</b> provides a way to check whether the variance
      between the samples is concomitatnt with their treatment groups<br> <br>

      <b>The Barchart</b> serves as an excellent way to
      summarize the up- and downregulated DEGs <br> <br>

      <b>Volcano</b> and <b>MA plots</b> are great for assessing DE results.
      A Volcano plot represents the relationship between significance and fold-change,
      whereas a MA plot shows the average expression of genes versus log fold-change <br> <br>

      <b>Heatmaps</b> can be used to visualize the top DEGs according to significance,
      assess expression patterns across the different conditions,
      and also as a quality control plot. <br> <br>
           </div>")
    } else{
      HTML("")
    }
  })



  observeEvent(input$figButton, {
    if (input$selectPlotCombo == "pca") {

      rv$plot <-
        plotPCA(rv$deTable[[2]], TRUE, rv$deTable[[3]], "treatment")


    } else if (input$selectPlotCombo == "scree") {
      rv$plot <- plotScree(rv$deTable[[2]])


    } else if (input$selectPlotCombo == "bar") {

      rv$barTable <- barTable(
        rv$merged,
        as.numeric(input$plotPvalue),
        as.numeric(input$plotFC)
      )

      rv$plot <- plotBarChart(rv$barTable)


    } else if (input$selectPlotCombo == "volcano") {
      rv$plot <- plotVP(
        rv$merged,
        as.numeric(input$plotFC),
        as.numeric(input$plotPvalue)
      )


    } else if (input$selectPlotCombo == "MA") {
      rv$plot <- plotMA(
        rv$merged,
        as.numeric(input$plotFC),
        as.numeric(input$plotPvalue)
      )


    } else if (input$selectPlotCombo == "heat") {
      rv$heatData <- plotHeatmapTop(
        rv$merged,
        as.numeric(input$topGeneNo),
        as.numeric(input$selectTypeHM),
        as.numeric(input$plotPvalue),
        as.numeric(input$plotFC),
        input$topNames,
        session
      )

      rv$plot <- rv$heatData[[1]]


    }

    output$mainPlot <- renderPlot({
      grid.draw(rv$plot)
    })
  })


  ## Download
  output$downloadTopHeat <- downloadHandler(
    filename = function() {
      paste(input$selectPlotCombo,
            "plot" ,
            device = ".png",
            sep = "")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(
          ...,
          width = width,
          height = height,
          units = "px",
          pointsize = 12
        )
      }
      ggsave(
        file,
        plot = rv$plot,
        device = device,
        width = 1280,
        height = 720,
        limitsize = FALSE
      )
    }
  )

  output$downloadTopGenes <- downloadHandler(
    filename = function() {
      paste("top", input$topGeneNo, "Genes" , ".csv", sep = "")
    },
    content = function(file) {
      data <- rv$heatData[[2]]

      write.csv(data, file)
    }
  )
}



#' Generate the corresponding table of the Barchart
#'
#' @param data Differential Expression results (deTable)
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @export
#' @return Returns the table to be used in the Barchart
barTable <- function(data, pvalue, fchange) {
  countTable <- data

  fchange <- log2(fchange)

  upSig.AB <-
    subset(countTable, FDR < pvalue & logFC > fchange)
  downSig.AB <-
    subset(countTable, FDR < pvalue & logFC < -fchange)

  upCount.AB <- nrow(upSig.AB)
  downCount.AB <- nrow(downSig.AB)

  comparison = c("", "")
  direction = c("Upregulated", "Downregulated")
  number_of_sig_genes = c(upCount.AB, downCount.AB)


  df <- data.frame(comparison,
                   direction,
                   number_of_sig_genes) # Merge vector into df

  return(df)
}

#' Generate the Barchart
#'
#' @param df The Table generated by the barTable Function
#' @export
#' @return Returns a Barchart
plotBarChart <- function(df) {
  p <-
    ggplot(data = df,
           aes(x = comparison, y = number_of_sig_genes, fill = direction)) +
    geom_bar(colour = "black",
             stat = "identity",
             position = "dodge") +
    ylab("Number of significant genes") +
    xlab("") +
    scale_fill_discrete(name = "Direction") +
    theme_classic(base_size = 13) +
    geom_text(
      aes(label = number_of_sig_genes),
      position = position_dodge(width = 0.9),
      vjust = 2
    )

  return(p)
}


#' Generate a Heatmap
#'
#' @param data Differential Expression results (deTable)
#' @param geneNo The number of genes to be displayed
#' @param type Filter by up-, down-, or absolute significe DE
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @param names Boolean - whether to show names or not
#' @param session Current R session
#'
#' @export
#' @return Returns a heatmap
plotHeatmapTop <-
  function(data,
           geneNo,
           type,
           pvalue,
           fchange,
           names,
           session) {

    out <- tryCatch(
      {

        fchange <- log2(fchange) # convert to log2


        #filter
        if (type == 1) {
          data <-
            subset(data, FDR < pvalue &
                     abs(data$logFC) > fchange) # absSig

        } else if (type == 2) {
          data <-
            subset(data, FDR < pvalue &
                     data$logFC > fchange) # upregSig

        } else if (type == 3) {
          data <-
            subset(data, FDR < pvalue &
                     data$logFC < -fchange) # downreg Sig

        } else if (type == 4) {
          data <-
            subset(data, FDR > pvalue |
                     abs(data$logFC) < fchange) # non-Sig

        }

        data <- data[order(data$FDR), ] # order by FDR
        data <- data[1:geneNo, ]
        data <- na.omit(data) # omit NANs


        # columns
        counts <- data[, 5:ncol(data)]   # Select  expression value columns
        counts.scaled <- t(scale(t(counts)))  # Convert FPKM to Z-scores

        row.order <-
          hclust(dist(counts.scaled), method = "average")$order # Cluster the data
        counts.scaled.clustered <-
          counts.scaled[row.order, ] # Order by row.order
        counts.scaled.clustered.m <-
          melt(as.matrix(counts.scaled.clustered)) # convert to ggplot-appropriate format

        hm.palette <-
          colorRampPalette(c("red", "white", "blue")) # Set the colour range


        plot <-
          ggplot(counts.scaled.clustered.m, aes(x = Var2, y = Var1, fill = value)) +
          geom_tile() +
          scale_fill_gradientn(colours = hm.palette(100), name ="Row Z-score") +
          ylab('Genes') + xlab('Samples') + theme_bw() +
          theme(
            axis.text.x = element_text(
              angle = 90,
              hjust = 1,
              size = 12
            ),
            axis.title = element_text(size = 14),
            axis.ticks = element_blank(),
            plot.background = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank()
          )

        if (!names) {
          plot <- plot + theme(axis.text.y = element_blank())
        }

        out <- list(plot, data)
      },
      error=function(cond) {
        sendSweetAlert(
          session = session,
          title = "No Genes to Plot",
          text = "Please lower parameter stringency",
          type = "error"
        )
        return()
      }
    )
    return(out)
  }


#' Generate a Volcano plot
#'
#' @param data Differential Expression results (deTable)
#' @param pValue P-value threshold
#' @param fcValue Fold-Change threshold
#' @export
#' @return Returns a Volcano plot
plotVP <- function(data, fcValue, pValue) {
  x <- na.omit(data)

  fcValue <- log2(fcValue)

  x$sig_flag <-
    as.factor(x$FDR < pValue & abs(x$logFC) > fcValue)

  VP <-
    ggplot(data = x, aes(x$logFC, y = -log10(x$FDR), colour = sig_flag)) +
    geom_point(size = 1.8) +
    xlab("-log10 Adjusted p-value") +
    ylab("Log2 Fold Change") +
    theme_classic() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    scale_colour_discrete(
      breaks = c("TRUE", "FALSE"),
      labels = c("Significant", "Non-significant")
    )

  return(VP)
}


#' Generate an MA plot
#'
#' @param data Differential Expression results (deTable)
#' @param pValue P-value threshold
#' @param fcValue Fold-Change threshold
#' @export
#' @return Returns a MA plot
plotMA <- function(data,
                   fcValue,
                   pValue) {

  x <- na.omit(data)

  fcValue <- log2(fcValue)

  exprValues <- x[, 5:ncol(data)]

  x$sig_flag <-
    as.factor(x$FDR < pValue & abs(x$logFC) > fcValue)
  x$mean_expression <- rowMeans(exprValues, na.rm = FALSE, dims = 1)


  ## Generate a MA plot
  MA <-
    ggplot(data = x, aes(
      x = log10(x$mean_expression + 0.001),
      y = x$logFC,
      colour = x$sig_flag
    )) +
    geom_point(size = 1.8) +
    geom_hline(aes(yintercept = 0), colour = "black", size = 0.75) +
    xlab("Log2 Mean Expression") +
    ylab("Log2 Fold Change") +
    theme_classic() +
    theme(legend.position = "bottom", legend.title = element_blank()) +
    scale_colour_discrete(
      breaks = c("TRUE", "FALSE"),
      labels = c("Significant", "Non-significant")
    )

  return(MA)
}
