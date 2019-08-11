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
          "Heatmap" = "heat",
          "Venn Diagram" = "venn"
        )
      ),



      #VP or MA
      conditionalPanel(
        condition = "input.selectPlotCombo == 'volcano' ||
                               input.selectPlotCombo == 'MA' ||
                               input.selectPlotCombo == 'bar' ||
                               input.selectPlotCombo == 'venn' ||
                               input.selectPlotCombo == 'heat'" ,
        ns = ns,

        numericInput(
          ns("plotPvalue"),
          label = ("FDR threshold"),
          value = 0.05,
          min = 0.0001,
          max = 0.5
        ),
        fluidRow(column(3, verbatimTextOutput(ns(
          "plotPvalue"
        )))),

        numericInput(
          ns("plotFC"),
          label = ("Fold change threshold"),
          value = 2,
          min = 0,
          max = 20
        ),

        verbatimTextOutput(ns("plotFC"))

      ),

      conditionalPanel(
        condition = "input.selectPlotCombo == 'volcano' ||
                               input.selectPlotCombo == 'MA'",
        ns = ns,

        selectInput(
          ns("selectConditionMA"),
          label = "Select DE Analysis Condition",
          choices = list(
            "A vs B" = 1,
            "B vs C" = 2,
            "A vs C" = 3
          ),
          selected = 1
        )


      ),


      # Heatmap
      conditionalPanel(
        condition = "input.selectPlotCombo == 'heat'",
        ns = ns,

        selectInput(
          ns("selectConditionTHM"),
          label = ("Select DE Analysis Condition"),
          choices = list(
            "All Conditions" = 0,
            "A vs B" = 1,
            "B vs C" = 2,
            "A vs C" = 3
          )
        ),

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

      ##### Download

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
        ns = ns,

        DT::dataTableOutput(ns("barTable"))
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
  plotChoices = list(
    "PCA plot" = "pca",
    "Scree plot" = "scree",
    "Barchart" = "bar",
    "Volcano Plot" = "volcano",
    "MA Plot" = "MA",
    "Heatmap" = "heat"
  )

  vennChoice = list("Venn Diagram" = "venn")


  output$helpPlotInfo <- renderUI({
    if(is.null(rv$plot)){
      HTML("<div style='border:2px solid blue; font-size: 14px; border-radius: 10px;'>

      <p style='text-align: center'><b> This tab enables the visualization of DE results. </b> </p> <br>

      <b>The Barchart</b> serves as an excellent way to summarize the up- and downregulated DEGs <br> <br>

      <b>Volcano</b> and <b>MA plots</b> are great ways to assess DE results.
      A Volcano plot represents the relationship between significance and fold-change,
      whereas a MA plot shows the average expression of genes versus log fold-change <br> <br>

      <b>Heatmaps</b> can be used to visualize the top DEGs according to significance,
      assess expression patterns across the different conditions,
      and also as a quality control plot. <br> <br>

      <b>Venn Diagrams</b> were implemented as an option to assess the overlapping DEGs between the comparisons when working with
      more than two groups.

           </div>")
    } else{
      HTML("")
    }
  })


  observeEvent(input$figButton, {
    if (input$selectPlotCombo == "pca") {
      print(head(rv$deTable))

      rv$plot <-
        plotPCA(rv$deTable, (rv$offset + 1):(ncol(rv$deTable)))


    } else if (input$selectPlotCombo == "scree") {
      rv$plot <- plotScree(rv$deTable, (rv$offset + 1):(ncol(rv$deTable)))


    } else if (input$selectPlotCombo == "bar") {
      rv$barTable <- barTable(
        rv$deTable,
        rv$conditionNo,
        as.numeric(input$plotPvalue),
        as.numeric(input$plotFC)
      )

      output$barTable <- DT::renderDataTable(DT::datatable(
        rv$barTable,
        colnames = c("Contrast", "Direction", "Significant Gene No")
      ))

      rv$plot <- plotBarChart(rv$barTable)


    } else if (input$selectPlotCombo == "volcano") {
      rv$plot <- plotVP(
        rv$deTable,
        as.numeric(input$selectConditionMA),
        as.numeric(input$plotFC),
        as.numeric(input$plotPvalue)
      )


    } else if (input$selectPlotCombo == "MA") {
      rv$plot <- plotMA(
        rv$deTable,
        rv$offset:((ncol(rv$deTable)) - 1),
        as.numeric(input$selectConditionMA),
        as.numeric(input$plotFC),
        as.numeric(input$plotPvalue)
      )


    } else if (input$selectPlotCombo == "heat") {
      rv$plot <- plotHeatmapTop(
        rv$deTable,
        rv$offset,
        rv$replicateNo,
        as.numeric(input$topGeneNo),
        as.numeric(input$selectTypeHM),
        as.numeric(input$plotPvalue),
        as.numeric(input$plotFC),
        as.numeric(input$selectConditionTHM),
        input$topNames
      )


      rv$heatData <- getHeatData(
        rv$deTable,
        rv$offset,
        rv$replicateNo,
        as.numeric(input$topGeneNo),
        as.numeric(input$selectTypeHM),
        as.numeric(input$plotPvalue),
        as.numeric(input$plotFC),
        as.numeric(input$selectConditionTHM)
      )


    } else if (input$selectPlotCombo == "venn") {
      rv$plot <-
        plotContrastVenn(rv$deTable, input$plotPvalue, input$plotFC)

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
      data <- rv$heatData

      write.csv(data, file)
    }
  )
}


#' Plot PCA
#'
#' Uses factoextra package to plot a PCA plot
#'
#' @param data Differential Expression results (deTable)
#' @param expressionColumns The columns with the Normalized Counts
#' @export
#' @return the PCA plot
plotPCA <- function(data, expressionColumns) {
  countTable <- data[, expressionColumns]
  countTable <- as.matrix(sapply(countTable, as.numeric))
  xx <- prcomp(t(countTable))

  return(fviz_pca_ind(xx, repel = TRUE))
}

#' Plot Scree
#'
#' Uses factoextra package to plot a Scree plot
#'
#' @param data Differential Expression results (deTable)
#' @param expressionColumns The columns with the Normalized Counts
#' @export
#' @return the Scree plot
plotScree <- function(data, expressionColumns) {
  x <- data[, expressionColumns]
  x <- as.matrix(sapply(x, as.numeric))
  xx <- prcomp(t(x))

  return(fviz_eig(xx))
}

#' Generate the corresponding table of the Barchart
#'
#' @param data Differential Expression results (deTable)
#' @param conditionNo The Fold-change column of the condition of interest
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @export
#' @return df The table to be used in the Barchart
barTable <- function(data, conditionNo, pvalue, fchange) {
  countTable <- data

  fchange <- log2(fchange)

  upSig.AB <-
    subset(countTable, FDR < pvalue & countTable[, 1] > fchange)
  downSig.AB <-
    subset(countTable, FDR < pvalue & countTable[, 1] < -fchange)

  upCount.AB <- nrow(upSig.AB)
  downCount.AB <- nrow(downSig.AB)

  if (conditionNo > 2) {
    # filter
    upSig.BC <-
      subset(countTable, FDR < pvalue & countTable[, 2] > fchange)
    upSig.AC <-
      subset(countTable, FDR < pvalue & countTable[, 3] > fchange)
    downSig.BC <-
      subset(countTable, FDR < pvalue & countTable[, 2] < -fchange)
    downSig.AC <-
      subset(countTable, FDR < pvalue & countTable[, 3] < -fchange)

    # gene No
    upCount.AC <- nrow(upSig.AC)
    upCount.BC <- nrow(upSig.BC)

    downCount.BC <- nrow(downSig.BC)
    downCount.AC <- nrow(downSig.AC)

  }


  if (conditionNo == 2) {
    comparison = c("A vs B", "A vs B")
    direction = c("Upregulated", "Downregulated")
    number_of_sig_genes = c(upCount.AB, downCount.AB)

  } else {
    comparison = c("A vs B", "A vs B", "B vs C", "B vs C", "A vs C", "A vs C")
    direction = c(
      "Upregulated",
      "Downregulated",
      "Upregulated",
      "Downregulated",
      "Upregulated",
      "Downregulated"
    )
    number_of_sig_genes = c(upCount.AB,
                            downCount.AB,
                            upCount.BC,
                            downCount.BC,
                            upCount.AC,
                            downCount.AC)
  }

  df <-
    data.frame(comparison, direction, number_of_sig_genes) # Merge vector into df

  # write.csv(df, file="output/DEGperCondition.csv", row.names = FALSE)

  return(df)
}

#' Generate the Barchart
#'
#' @param df The Table generated by the barTable Function
#' @export
#' @return the Barchart
plotBarChart <- function(df) {
  p <-
    ggplot(data = df,
           aes(x = comparison, y = number_of_sig_genes, fill = direction)) +
    geom_bar(colour = "black",
             stat = "identity",
             position = "dodge") +
    ylab("Number of significant genes") +
    xlab("Contrast") +
    scale_fill_discrete(name = "Direction") +
    theme_classic(base_size = 13) +
    geom_text(
      aes(label = number_of_sig_genes),
      position = position_dodge(width = 0.9),
      vjust = 2
    )

  return(p)
}

#' Generate a Venn Diagram
#'
#' Creates a Venn diagram of the DE genes between the different conditions
#' @param data Differential Expression results (deTable)
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @export
#' @return The Venn Diagram
plotContrastVenn <- function(data, pvalue, fchange) {
  countTable <- data

  fchange <- log2(fchange)

  x1_sig <-
    subset(countTable, FDR < pvalue & abs(countTable[, 1]) > fchange)
  x2_sig <-
    subset(countTable, FDR < pvalue & abs(countTable[, 2]) > fchange)
  x3_sig <-
    subset(countTable, FDR < pvalue & abs(countTable[, 3]) > fchange)


  AvB <- x1_sig[1]
  BvC <- x2_sig[1]
  AvC <- x3_sig[1]

  xlist <- c(AvB, BvC, AvC)

  names(xlist) <- c("B vs A", "B vs C", "C vs A")

  venn.plot <- venn.diagram(
    xlist,
    filename = NULL,
    fill = c("red", "blue", "green"),
    lty = "blank",
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.default.pos = "outer"
  )
  grid.newpage()
  grid.draw(venn.plot)

  return(venn.plot)
}

#' Generate a Heatmap
#'
#' @param data Differential Expression results (deTable)
#' @param offset The last column of DE results
#' @param replicates The number of replicates
#' @param geneNo The number of genes to be displayed
#' @param type Filter by up-, down-, or absolute significe DE
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @param conditionFC The Fold-change column of the condition of interest
#' @param names Boolean - whether to show names or not
#' @export
#' @return The heatmap
plotHeatmapTop <-
  function(data,
           offset,
           replicates,
           geneNo,
           type,
           pvalue,
           fchange,
           conditionFC,
           names) {
    table <- data

    fchange <- log2(fchange) # convert to log2

    # save appropriate columns
    if (conditionFC == 1) {
      columns = (offset + 1):((offset) + 2 * replicates) #AvB
    } else if (conditionFC == 2) {
      columns = ((offset + 1) + replicates):((offset) + 3 * replicates) #BvC
    } else if (conditionFC == 3) {
      columns = c((offset + 1):(offset + replicates),
                  ((offset + 1) + 2 * replicates):((offset) + 3 * replicates)) #AvC
    } else{
      conditionFC = c(1, 2, 3)
      columns = ((offset) + 1):(ncol(table))
    }

    #filter
    if (type == 1) {
      table <-
        subset(table, FDR < pvalue &
                 abs(table[, conditionFC]) > fchange) # absSig

    } else if (type == 2) {
      table <-
        subset(table, FDR < pvalue &
                 table[, conditionFC] > fchange) # upregSig

    } else if (type == 3) {
      table <-
        subset(table, FDR < pvalue &
                 table[, conditionFC] < -fchange) # downreg Sig

    } else if (type == 4) {
      table <-
        subset(table, FDR > pvalue |
                 abs(table[, conditionFC]) < fchange) # non-Sig

    }

    table <- table[order(table$FDR), ] # order by FDR
    table <- table[1:geneNo, ]

    table <- na.omit(table) # omit NANs

    counts <- table[, columns]   # Select  expression value columns
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
      geom_tile() +  scale_fill_gradientn(colours = hm.palette(100), name =
                                            "Row Z-score") +
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

    return(plot)
  }


#' Returns Heatmap data
#'
#' @param data Differential Expression results (deTable)
#' @param offset The last column of DE results
#' @param replicates The number of replicates
#' @param geneNo The number of genes to be displayed
#' @param type Filter by up-, down-, or absolute significe DE
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @param conditionFC The Fold-change column of the condition of interest
#' @param names Boolean - whether to show names or not
#' @export
#' @return table The data that corresponds to the displayed heatmap
getHeatData <-
  function(data,
           offset,
           replicates,
           geneNo,
           type,
           pvalue,
           fchange,
           conditionFC) {
    table <- data

    fchange <- log2(fchange) # convert to log2

    # save appropriate columns
    if (conditionFC == 1) {
      columns = (offset + 1):((offset) + 2 * replicates) #AvB
    } else if (conditionFC == 2) {
      columns = ((offset + 1) + replicates):((offset) + 3 * replicates) #BvC
    } else if (conditionFC == 3) {
      columns = c((offset + 1):(offset + replicates),
                  ((offset + 1) + 2 * replicates):((offset) + 3 * replicates)) #AvC
    } else{
      conditionFC = c(1, 2, 3)
      columns = ((offset) + 1):(ncol(table))
    }




    #filter
    if (type == 1) {
      table <-
        subset(table, FDR < pvalue &
                 abs(table[, conditionFC]) > fchange) # absSig

    } else if (type == 2) {
      table <-
        subset(table, FDR < pvalue &
                 table[, conditionFC] > fchange) # upregSig

    } else if (type == 3) {
      table <-
        subset(table, FDR < pvalue &
                 table[, conditionFC] < -fchange) # downreg Sig

    } else if (type == 4) {
      table <-
        subset(table, FDR > pvalue |
                 abs(table[, conditionFC]) < fchange) # non-Sig

    }

    table <- table[order(table$FDR), ] # order by FDR
    table <- table[1:geneNo, ]

    table <- na.omit(table) # omit NANs

    return(table)
  }


#' Generate a Volcano plot
#'
#' @param data Differential Expression results (deTable)
#' @param pValue P-value threshold
#' @param fcValue Fold-Change threshold
#' @param fcColumn The Fold-change column of the condition of interest
#' @export
#' @return The Volcano plot
plotVP <- function(data, fcColumn, fcValue, pValue) {
  x <- data

  fcValue <- log2(fcValue)

  x$sig_flag <-
    as.factor(x$FDR < pValue & abs(x[, fcColumn]) > fcValue)

  x <- na.omit(x)

  VP <-
    ggplot(data = x, aes(x[, fcColumn], y = -log10(x$FDR), colour = sig_flag)) + geom_point(size =
                                                                                              1.8) +
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
#' @param expColumns normalized count columns of condition of interest
#' @param pValue P-value threshold
#' @param fcValue Fold-Change threshold
#' @param fcColumn The Fold-change column of the condition of interest
#' @export
#' @return The MA plot
plotMA <- function(data,
                   expColumns,
                   fcColumn,
                   fcValue,
                   pValue) {
  x <- data

  fcValue <- log2(fcValue)

  x <- na.omit(x)

  exprValues <- x[, expColumns]

  x$sig_flag <-
    as.factor(x$FDR < pValue & abs(x[, fcColumn]) > fcValue)
  x$mean_expression <- rowMeans(exprValues, na.rm = FALSE, dims = 1)


  ## Generate a MA plot
  MA <-
    ggplot(data = x, aes(
      x = log10(x$mean_expression + 0.001),
      y = x[, fcColumn],
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
