#' Bulk Functional Annotation Tab UI
#'
#' @export
#' @return None
bulk_goDataUI <- function(id) {
  ns <- NS(id)

  tagList(# Sidebar panel for inputs ----
          sidebarPanel(tabsetPanel(
            id = ns("goSideTabSet"),
            tabPanel(
              value = ("getGOGenesTab"),
              title = "Get DE Genes",

              numericInput(
                ns("goGetPvalue"),
                label = "Corrected p-value threshold",
                value = 0.05,
                min = 0.00001,
                max = 0.5
              ),
              fluidRow(column(3, verbatimTextOutput("goGetPvalue"))),


              numericInput(
                ns("goGetFC"),
                label = "Fold change threshold",
                value = 2,
                min = 0,
                max = 20
              ),
              fluidRow(column(3, verbatimTextOutput("goGetFC"))),


              radioButtons(
                ns("goGetType"),
                label = "Direction",
                c(
                  "Differentially Expressed" = 1,
                  "Upregulated" = 2,
                  "Downregulated" = 3
                )
              ),


              radioButtons(
                ns("goGetCondition"),
                label = "Condition",
                c(
                  "All Conditions" = 0 ,
                  "A vs B" = 1,
                  "B vs C" = 2,
                  "A vs C" = 3
                )
              ),


              # Input: Button gets Differentially Expressed Genes ----
              actionButton(ns("goGetButton"), label = "Get DE genes")

            ),

            tabPanel(
              value = ("getGOTermsTab"),
              title = "Get GO Terms",

              numericInput(
                ns("goTermPvalue"),
                label = "P-value threshold",
                value = 0.05,
                min = 0.00001,
                max = 0.5
              ),
              fluidRow(column(3, verbatimTextOutput(
                ns("goTermPvalue")
              ))),

              checkboxInput(ns("goTermFDR"), label = "FDR Correction", value = TRUE),

              radioButtons(
                ns("goTermSymbol"),
                label = "Symbol Type",
                c("Gene Symbol" = "geneSymbol", "Ensembl gene ID" = "ensGene")
              ),

              radioButtons(
                ns("goTermTest"),
                label = "Test Type",
                c(
                  "Cellular Component" = "GO:CC",
                  "Biological Process" = "GO:BP",
                  "Molecular Function" = "GO:MF",
                  "KEGG Pathways" = "KEGG"
                )
              ),


              selectInput(
                ns("goTermGenome"),
                label = "Select Genome",
                choices = list(
                  "Homo Sapiens (hg19)" = "hg19",
                  "Homo Sapiens (hg18)" = "hg18",
                  "Homo Sapiens (hg17)" = "hg17",
                  "Mus musculus (mm9)" = "mm9",
                  "Mus musculus (mm9)" = "mm8",
                  "Mus musculus (mm7)" = "mm7"
                )
              ),

              actionButton(ns("goTermButton"), label = "Get GO Terms"),
              actionButton(ns("goHistButton"), label = "Plot Histogram"),

              tags$hr(),

              textInput(ns("goInfoInput"), "Enter GO ID"),
              actionButton(ns("goInfoButton"), label = "Get Information")

            )
          )),

          # Main panel for displaying outputs ----
          mainPanel(
            tabsetPanel(
              id = ns("goMainTabSet"),
              tabPanel(
                value = ("goGenesTab"),
                title = "DE Genes Table",
                verbatimTextOutput(ns("goGenesText"), placeholder = T),

                tags$hr(),

                DT::dataTableOutput(ns("goGenesTable"))
              ),

              tabPanel(
                value = ("goTableTab"),
                title = "GO Term Table",
                DT::dataTableOutput(ns("goTermTable"))
              ),

              tabPanel(
                value = ("goHistTab"),
                title = "GO Term Histogram",
                checkboxInput(ns("goHistCheck"), "Function Names", FALSE),
                plotOutput(ns("goHistPlot"))
              ),

              tabPanel(
                value = ("goInfoTab"),
                title = "GO Term Info",
                verbatimTextOutput(ns("goInfoText"))
              )
            )
          ))
}


#' Bulk Functional Annotation Tab Server
#'
#' @param counts Unfiltered Count Table (Reactive Value)
#' @param de Differential Expression Results (Reactive Value)
#' @return None
bulk_goData <- function(input, output, session, counts, de) {
  rv <- reactiveValues()


  # Functional Annotation ------
  observeEvent(input$goGetButton, {
    rv$goGetGenes <-
      getDEgenes(
        de$deTable,
        input$goGetType,
        input$goGetPvalue,
        input$goGetFC,
        as.numeric(input$goGetCondition)
      )

    x <- as.data.frame(rv$goGetGenes)

    output$goGenesTable <-
      DT::renderDataTable(x, colnames = ("Differentially Expressed Genes"))

    output$goGenesText <-
      renderText({
        paste("Number of DE genes:", nrow(as.data.frame(x)), sep = " ")
      })

    updateTabsetPanel(session, "goSideTabSet",
                      selected = "getGOTermsTab")

  })

  observeEvent(input$goTermButton, {
    rv$goTermTable <-
      runGOSEQ(
        rv$goGetGenes,
        counts$countTable,
        input$goTermGenome,
        input$goTermSymbol,
        input$goTermTest,
        input$goTermFDR,
        input$goTermPvalue,
        2, # Bulk app
        session
      )


    output$goTermTable <- DT::renderDataTable(DT::datatable(
      rv$goTermTable,
      rownames = FALSE,
      colnames = c(
        "ID",
        "Overrepresented Genes p-value",
        "Overrepresented Genes p-value",
        "DE Genes No",
        "Total Gene No",
        "GO Term",
        "Ontology"
      )
    ))

    updateTabsetPanel(session, "goMainTabSet", selected = "goTableTab")

  })



  observeEvent(input$goHistButton, {
    rv$hist <- histGoTerms(rv$goTermTable, input$goHistCheck)

    output$goHistPlot <- renderPlot({
      rv$hist

    })

    # ggsave("figures/GOTermHistogram.png", plot = rv$hist, device = png(),
    # width = 12, height = 8, limitsize = FALSE)

    updateTabsetPanel(session, "goMainTabSet", selected = "goHistTab")
  })

  observeEvent(input$goInfoButton, {
    go.term <- GOTERM[[input$goInfoInput]]

    output$goInfoText <- renderText({
      paste(GOID(go.term),
            Term(go.term),
            Definition(go.term),
            Synonym(go.term),
            sep = "\n")
    })

    updateTabsetPanel(session, "goMainTabSet", selected = "goInfoTab")
  })
}


#' Get Differentially Expressed Genes (Bulk)
#'
#' Filters the DE results table and returns the names of the DE genes
#'
#' @param data Differential Expression Reslts Table
#' @param type Filters by abs. significant, upregulation or downregulation
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @param condition Filters by a given condition (AvB, BvC, etc.)
#' @return gene.vector a vector with DE gene names
getDEgenes <- function(data, type, pvalue, fchange, condition) {
  table <- data

  fchange <- log2(fchange) # convert to log2


  if (condition == 0) {
    #to filter for all conditions
    condition = c(1, 2, 3)
  }


  #filter
  if (type == 1) {
    table <-
      subset(table, FDR < pvalue &
               abs(table[, condition]) > fchange) # absSig

  } else if (type == 2) {
    table <-
      subset(table, FDR < pvalue &
               table[, condition] > fchange) # upregSig

  } else if (type == 3) {
    table <-
      subset(table, FDR < pvalue &
               table[, condition] < -fchange) # downreg Sig

  }

  table <- na.omit(table) # omit NANs

  gene.vector <- as.vector(row.names(table))

  # write.csv(gene.vector, file="output/GO_DEgenes.csv", row.names = FALSE)

  return(gene.vector)
}

#' Functional Annotation Pipeline (Bulk)
#'
#' Runs the GO pipeline according to the users preferences
#'
#' @param deGenes a vector with DE gene names
#' @param data Filters by abs. significant, upregulation or downregulation
#' @param genome The genome to be used
#' @param symbol Type of Gene symbol
#' @param testType Type of Test (Cellular component, Molecular function, etc.)
#' @param pCorrect Whether to filter by FDR p-value or uncorrected p-value
#' @param pValue P-value threshold
#' @param appNum Format entered data appropriately to the app launched
#'
#' @return GO.wall A dataframe with Functional Annotation Results
runGOSEQ <-
  runGOSEQ <-
  function(deGenes,
           data,
           genome,
           symbol,
           testType,
           pCorrect,
           pValue,
           appNum,
           session) {

    out <- tryCatch(
      {
        if(appNum == 1){ #Single-cell
          assayed.genes <- as.vector(row.names(data)) #* countTable to data
        } else if(appNum == 2){ # Bulk
          countTable  <- data[, -1]
          rownames(countTable) <- data[, 1]
          assayed.genes <- as.vector(row.names(countTable))
        }

        gene.vector = as.integer(assayed.genes %in% deGenes)
        names(gene.vector) = assayed.genes

        pwf = nullp(gene.vector, genome, symbol)

        GO.wall = goseq(pwf, genome, symbol, test.cats = testType)

        if (pCorrect) {
          GO.wall$over_represented_pvalue <-
            p.adjust(GO.wall$over_represented_pvalue, method = "BH")
        }

        GO.wall <-
          subset(GO.wall, GO.wall$over_represented_pvalue < pValue)

        if (testType == "GO:CC") {
          holder = "CellularComponent"
        } else if (testType == "GO:BP") {
          holder = "BiologicalFunction"
        } else if (testType == "GO:MF") {
          holder = "MolecularFunction"
        } else{
          holder = "KEGG"
        }

        # write.csv(GO.wall, file = paste("output/GO_Terms_", holder, ".csv", sep=""), row.names = FALSE)
        return(GO.wall)
      },
      error=function(cond) {

        sendSweetAlert(
          session = session,
          title = "Incorrect Gene type",
          text = "Ensure that the correct Gene symbol or genome is chosen",
          type = "error"
        )

        return()      # Choose a return value in case of error
      }
    )
  }

#' Generate GO Histogram
#'
#' Creates a GO Histogram with GO terms or GO IDs
#'
#' @param data Functional Annotation results
#' @param gof Boolean to use GO:Terms or GO:ID
#' @return p The histogram
histGoTerms <- function(data, gof) {
  if (nrow(data) > 10) {
    data <- head(data, n = 10)
  }

  if (gof) {
    goVector <- as.list(NULL)
    n <- length(data$category)

    for (i in 1:n) {
      go.term <- (GOTERM[data$category[i]])
      goVector[i] <- Term(go.term)

    }

    data$goVector <- goVector
    out <- transform(data, goVector = unlist(goVector))
    data$category <- as.factor(out$goVector)
    levels(data$category) <-
      gsub("(.{20,}?)\\s", "\\1\n", levels(data$category))
  }


  p <- ggplot(data, aes(category,-log2(over_represented_pvalue))) +
    geom_col(color = "black") +
    theme_bw() +
    labs(y = "-log2(p-value)", x = "") +
    theme(
      text = element_text(size = 20),
      axis.text.x = element_text(
        angle = 90,
        hjust = 0,
        vjust = 0.2
      ),
      axis.ticks.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_blank()
    )

  # png("figures/GOTermHistogram.png", height=900, width=1200)
  print(p)
  dev.off()

  return(p)
}
