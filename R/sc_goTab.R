#' SC Functional Annotation Tab UI
#'
#' @export
#' @return None
sc_goUI <- function(id) {
  ns <- NS(id)

  tagList(# Sidebar panel for inputs ----
          sidebarPanel(tabsetPanel(
            id = ns("goSideTabSet"),
            tabPanel(
              value = "getGOGenesTab",
              title = "Get DE Genes",

              checkboxInput(ns("goClustCheck"), label = h4("Get Genes from All Clusters"), FALSE),

              conditionalPanel(
                condition = "!input.goClustCheck",
                ns = ns,

                uiOutput("goClustNoInputUI")
              ),

              numericInput(
                ns("goGetPvalue"),
                label = "Adjusted P-value Threshold <",
                value = 0.05,
                min = 0.00001,
                max = 0.5
              ),

              fluidRow(column(3, verbatimTextOutput(
                ns("goGetPvalue")
              ))),


              numericInput(
                ns("goGetFC"),
                label = "Fold-Change Threshold >",
                value = 2,
                min = 0,
                max = 20
              ),

              fluidRow(column(3, verbatimTextOutput(ns(
                "goGetFC"
              )))),


              radioButtons(
                ns("goGetType"),
                label = "Direction",
                c(
                  "Differentially Expressed" = 1,
                  "Upregulated" = 2,
                  "Downregulated" = 3
                )
              ),

              actionButton(ns("goGetButton"), label = "Get DE genes")
            ),

            tabPanel(
              value = "getGOTermsTab",
              title = "Get GO Terms",

              numericInput(
                ns("goTermPvalue"),
                label = "P-value threshold <",
                value = 0.05,
                min = 0.00001,
                max = 0.5
              ),
              fluidRow(column(3, verbatimTextOutput(
                ns("goTermPvalue")
              ))),

              checkboxInput(ns("goTermFDR"),
                            label = "FDR Correction",
                            value = TRUE),

              radioButtons(
                ns("goTermSymbol"),
                label = "Symbol Type",
                c("Gene Symbol" = "geneSymbol", "Ensembl gene ID" = "ensGene")
              ),

              radioButtons(
                ns("goTermTest"),
                label = "Ontology Type",
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
                  "Homo sapiens (hg19)" = "hg19",
                  "Homo sapiens (hg18)" = "hg18",
                  "Homo sapiens (hg17)" = "hg17",
                  "Mus musculus (mm9)" = "mm9",
                  "Mus musculus (mm9)" = "mm8",
                  "Mus musculus (mm7)" = "mm7",
                  "Danio rerio (danRer5)" = "danRer5",
                  "Drosophila melanogaster (dm3)" = "dm3",
                  "E. coli K12" = "E. coli K12"
                )
              ),

              actionButton(ns("goTermButton"), label = "Get GO Terms"),

              tags$hr(),

              tags$b("Top 10 GO Term Histogram"),
              checkboxInput(ns("goHistCheck"), "Function Names", TRUE),
              actionButton(ns("goHistButton"), label = "Plot Histogram"),

              tags$hr(),

              textInput(ns("goInfoInput"), "Enter GO:ID of interest:"),
              actionButton(ns("goInfoButton"), label = "Get Information")

            )
          )),

          # Main panel for displaying outputs ----
          mainPanel(
            tabsetPanel(
              id = ns("goMainTabSet"),
              tabPanel(
                value = "goGenesTab",
                title = "DE Genes Table",

                htmlOutput(ns("helpGoGeneInfo")),

                verbatimTextOutput(ns("goGenesText"), placeholder = T),

                tags$hr(),

                DT::dataTableOutput(ns("goGenesTable"))
              ),

              tabPanel(
                value = "goTableTab",
                title = "GO Term Table",
                DT::dataTableOutput(ns("goTermTable")),
                conditionalPanel(condition = "input.goTermButton > 0",
                                 ns = ns,
                                 downloadButton(ns("goDownload"), "Download Table")
                )
              ),

              tabPanel(
                value = "goHistTab",
                title = "GO Term Histogram",
                plotOutput(ns("goHistPlot"), width = "800px", height = "500px")
              ),

              tabPanel(
                value = "goInfoTab",
                title = "GO Term Info",
                verbatimTextOutput(ns("goInfoText"))
              )
            )
          ))
}


#' SC Functional Annotation Tab Server
#'
#' @param countsT Unfiltered Count Table (Reactive Value)
#' @param de Differential Expression Results (Reactive Value)
#' @return None
sc_go <- function(input, output, session, de, countsT) {
  go <- reactiveValues()

  output$helpGoGeneInfo <- renderUI({
    if(input$goGetButton == 0){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px;'>
        The <i> 'Get DE Genes' </i> tab enables DE genes
        obtained from the DE analysis to be pre-filtered according to: <br>
        Fold-change, adj. P-value threshold,
        and according to a specific cluster <br>
        Once the DE genes are filtered, the 'Get GO Terms'
           tab will be automatically selected. </div>")
    } else {
      if(is.null(go$goTermTable)){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bot: 8px; border-radius: 10px;'>
      The <i> 'Get GO Terms' </i> tab provides a comprehensive
      Functional Annotation Pipeline using the filtered DE genes. <br> <br>

      Prior to running the pipeline, please specify the following parameters:
      adjusted or non-adjusted P-value threshold for ontologies;
      Symbol type - the type of symbols used for the genes in the count table;
      Ontology type; and Genome of interest. <br> <br>

      Once a table with results is returned,
      proceed to visualizing the top 10 GO terms results via the histogram option
      and to exploring GO terms of interest using their GO:IDs </div>" )
      }else{
        HTML("")
      }

    }
  })

  observeEvent(input$goGetButton, {
    if(!is.null(de$markers)){
      go$goGetGenes <-
        sc_getDEgenes(
          de$markers,
          input$goGetType,
          input$goGetPvalue,
          input$goGetFC,
          input$goClustCheck,
          input$goClustNoInput
        )

      x <- as.data.frame(go$goGetGenes)

      output$goGenesTable <-
        DT::renderDataTable(x, colnames = ("Differentially Expressed Genes"))

      output$goGenesText <-
        renderText({
          paste("Number of DE genes:", nrow(as.data.frame(x)), sep = " ")
        })

      updateTabsetPanel(session,
                        "goSideTabSet",
                        selected = "getGOTermsTab")
    }
  })



  observeEvent(input$goTermButton, {

    if(!is.null(go$goGetGenes)){

      waiter_show(tagList(spin_folding_cube(), h2("Loading ...")))

      go$goTermTable <-
        runGOSEQ(
          go$goGetGenes,
          countsT$countTable,
          input$goTermGenome,
          input$goTermSymbol,
          input$goTermTest,
          input$goTermFDR,
          input$goTermPvalue,
          1, # Single-Cell
          session
        )

      output$goTermTable <- DT::renderDataTable(
        go$goTermTable %>% datatable() %>%
          formatSignif(columns = c(2:3), digits = 4),
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
      )

      waiter_hide()

      updateTabsetPanel(session, "goMainTabSet", selected = "goTableTab")
    }
  })


  observeEvent(input$goHistButton, {
    if(!is.null(go$goTermTable)){
      go$goHistPlot <- histGoTerms(go$goTermTable, input$goHistCheck, session)

      output$goHistPlot <- renderPlot({
        go$goHistPlot
      })


      updateTabsetPanel(session, "goMainTabSet", selected = "goHistTab")
    }
  })

  observeEvent(input$goInfoButton, {

    output$goInfoText <- renderText({
      goInfo(input$goInfoInput, session)
    })

    updateTabsetPanel(session, "goMainTabSet", selected = "goInfoTab")
  })

  output$goDownload <- downloadHandler(
    filename = function() {
      paste(format(Sys.time(), "%y-%m-%d_%H-%M"), "_goTermResults" , ".csv", sep = "")
    },
    content = function(file) {
      data <- go$goTermTable

      write.csv(data, file)
    }
  )
}




#' Get Differentially Expressed Genes (Single Cell)
#'
#' Filters the DE results table and returns the names of the DE genes
#'
#' @param data Differential Expression Reslts Table
#' @param type Filters by abs. significant, upregulation or downregulation
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @param condition Filters by a given condition (AvB, BvC, etc.)
#' @return Returns a vector with DE gene names
sc_getDEgenes <-
  function(data,
           type,
           pvalue,
           fchange,
           byClust,
           clusterNo) {
    table <- data
    colnames(table)[5] <- "FDR"

    fchange <- log(fchange)


    if (!byClust) {
      # look only at a particular cluster
      table <- table[table$cluster == clusterNo, ]

    }

    #filter
    if (type == 1) {
      table <-
        subset(table, FDR < pvalue & abs(avg_logFC) > fchange) # absSig

    } else if (type == 2) {
      table <- subset(table, FDR < pvalue &
                        avg_logFC > fchange) # upregSig

    } else if (type == 3) {
      table <-
        subset(table, FDR < pvalue & avg_logFC < -fchange) # downreg Sig
    }

    table <- na.omit(table) # omit NANs

    gene.vector <- as.vector(row.names(table))

    return(gene.vector)
  }
