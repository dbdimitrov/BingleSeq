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

              h4("Filter DE Genes"),

              numericInput(
                ns("goGetPvalue"),
                label = "Adjusted P-value threshold <",
                value = 0.05,
                min = 0.00001,
                max = 0.5
              ),
              fluidRow(column(3, verbatimTextOutput("goGetPvalue"))),


              numericInput(
                ns("goGetFC"),
                label = "Fold-change threshold >",
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

              actionButton(ns("goGetButton"), label = "Get DE genes")

            ),

            tabPanel(
              value = ("getGOTermsTab"),
              title = "Get GO Terms",

              h4("Get Gene Ontologies"),


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
                value = ("goGenesTab"),
                title = "DE Genes Table",

                htmlOutput(ns("helpGoGeneInfo")),

                verbatimTextOutput(ns("goGenesText"), placeholder = T),

                tags$hr(),

                DT::dataTableOutput(ns("goGenesTable"))
              ),

              tabPanel(
                value = ("goTableTab"),
                title = "GO Term Table",
                DT::dataTableOutput(ns("goTermTable")),
                conditionalPanel(condition = "input.goTermButton > 0",
                                 ns = ns,
                                 downloadButton(ns("goDownload"), "Download Table")
                )
              ),

              tabPanel(
                value = ("goHistTab"),
                title = "GO Term Histogram",
                plotOutput(ns("goHistPlot"), width = "800px", height = "700px")
              ),

              tabPanel(
                value = ("goInfoTab"),
                title = "GO Term Info",
                verbatimTextOutput(ns("goInfoText"))
              )
            )
          ))
}



#' Bulk Gene Ontology Tab Server
#'
#' @param counts Unfiltered Count Table (Reactive Value)
#' @param de Differential Expression Results (Reactive Value)
#' @return None
bulk_goData <- function(input, output, session, counts, de) {
  go <- reactiveValues()

  output$helpGoGeneInfo <- renderUI({
    if(input$goGetButton == 0){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px'>
        
        
        This tab enables the Gene Ontology and Pathway analysis of deregulated 
        gene sets using the over representation analysis approach implemented 
        within the GOseq package.
        <br> <br>
        
        The <i> 'Get DE Genes' </i> subtab enables
        DE genes obtained from the DE analysis to be pre-filtered according to:
        Fold-change, adj. P-value threshold <br>
        Once the DE genes are filtered,
           the 'Get GO Terms' tab will be automatically selected. </div>")
    } else {
      if(is.null(go$goTermTable)){
        HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px;'>
        The <i> 'Get GO Terms' </i> subtab provides a comprehensive
        Functional Annotation Pipeline using the filtered DE genes. <br> <br>

        Prior to running the pipeline, please specify the following parameters:
        adjusted or non-adjusted P-value threshold for ontologies;
        Symbol type - the type of symbols used for the genes in the count table;
        Ontology type; and Genome of interest. <br> <br>

        Once a table with results is returned, proceed to visualizing
        the top 10 GO terms results via the histogram option
        and to exploring GO terms of interest using their GO:IDs </div>" )
      }else{
        HTML("")
      }

    }
  })

  # Functional Annotation ------
  observeEvent(input$goGetButton, {

    go$goGetGenes <-
      getDEgenes(
        de$merged,
        input$goGetType,
        input$goGetPvalue,
        input$goGetFC
      )

    output$goGenesTable <-
      DT::renderDataTable(as.data.frame(go$goGetGenes),
                          colnames = ("Differentially Expressed Genes"))

    output$goGenesText <-
      renderText({
        paste("Number of DE genes:", nrow(as.data.frame(go$goGetGenes)), sep = " ")
      })

    updateTabsetPanel(session, "goSideTabSet",
                      selected = "getGOTermsTab")

  })

  observeEvent(input$goTermButton, {

    waiter_show(html=tagList(spin_folding_cube(), h2("Loading ...")))

    go$goTermTable <-
      runGOSEQ(
        go$goGetGenes,
        counts$countTable,
        input$goTermGenome,
        input$goTermSymbol,
        input$goTermTest,
        input$goTermFDR,
        input$goTermPvalue,
        2, # Bulk app
        session
      )


    output$goTermTable <- DT::renderDataTable(
      go$goTermTable %>% datatable(),
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

  })

  observeEvent(input$goHistButton, {
    go$hist <- histGoTerms(go$goTermTable, input$goHistCheck, session)

    output$goHistPlot <- renderPlot({
      go$hist
    })

    updateTabsetPanel(session, "goMainTabSet", selected = "goHistTab")
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


#' Get Differentially Expressed Genes (Bulk)
#'
#' Filters the DE results table and returns the names of the DE genes
#'
#' @param data Differential Expression Reslts Table
#' @param type Filters by abs. significant, upregulation or downregulation
#' @param pvalue P-value threshold
#' @param fchange Fold-Change threshold
#' @return Returns a vector with DE gene names
getDEgenes <- function(data, type, pvalue, fchange) {
  
  fchange <- log2(fchange) # convert to log2
  
  #filter
  if (type == 1) {
    table <-
      subset(data, FDR < pvalue &
               abs(logFC) > fchange) # absSig
    
  } else if (type == 2) {
    table <-
      subset(data, FDR < pvalue &
               logFC > fchange) # upregSig
    
  } else if (type == 3) {
    table <-
      subset(data, FDR < pvalue &
               logFC < -fchange) # downreg Sig
    
  }
  
  table <- na.omit(table) # omit NANs
  
  gene.vector <- as.vector(row.names(table))
  
  return(gene.vector)
}



#' Get GO gene sets for GOSeq
#'
#' Filters the DE results table and returns the names of the DE genes
#' @param assayed.genes Genes in the DE set
#' @param genome Model organism
#' @param geneset_symbol Gene nomenclature 
#' 
#' @return Returns a data frame with two columns: genes and GO
getGenesetsGO <- function(assayed.genes, genome, geneset_symbol){
  if (startsWith(genome, "hg")){
    require(org.Hs.eg.db)
    require(AnnotationDbi)
    go_data <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys=assayed.genes,
                                     keytype = geneset_symbol, 
                                     columns = c(geneset_symbol, "GO")) %>%
      dplyr::select(all_of(geneset_symbol), "GO")
    
  } else if(startsWith(genome, "mm")) {
    require(org.Mm.eg.db)
    go_data <- AnnotationDbi::select(org.Hs.eg.db,
                                     keys=assayed.genes,
                                     keytype = geneset_symbol, 
                                     columns = c(geneset_symbol, "GO")) %>%
      dplyr::select(all_of(geneset_symbol), "GO")
    
  } else if(startsWith(genome, "dan")) {
    require(org.Dr.eg.db)
    go_data <- AnnotationDbi::select(org.Dr.eg.db,
                                     keys=assayed.genes,
                                     keytype = geneset_symbol, 
                                     columns = c(geneset_symbol, "GO")) %>%
      dplyr::select(all_of(geneset_symbol), "GO")
    
  } else if(startsWith(genome, "dm")) {
    require(org.Dr.eg.db)
    go_data <- AnnotationDbi::select(org.Dm.eg.db,
                                     keys=assayed.genes,
                                     keytype = geneset_symbol, 
                                     columns = c(geneset_symbol, "GO")) %>%
      dplyr::select(all_of(geneset_symbol), "GO")
   } else if(startsWith(genome, "E.")) {
      require(org.Dr.eg.db)
      go_data <- AnnotationDbi::select(org.EcK12.eg.db,
                                       keys=assayed.genes,
                                       keytype = geneset_symbol, 
                                       columns = c(geneset_symbol, "GO")) %>%
        dplyr::select(all_of(geneset_symbol), "GO")
    } 
  
  return(go_data)
}




#' Functional Annotation Pipeline (Both pipelines)
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
#' @return Returns a dataframe with Functional Annotation Results
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
        
        geneset_symbol = switch(symbol,
                                "ensGene" = "ENSEMBL",
                                "geneSymbol" = "SYMBOL")
        
        go_data <- getGenesetsGO(assayed.genes, genome, geneset_symbol)
        
        pwf = nullp(gene.vector, genome, symbol)
        
        
        GO.wall = goseq(pwf, gene2cat = go_data, test.cats = testType)
        
        if (pCorrect) {
          GO.wall$over_represented_pvalue <-
            p.adjust(GO.wall$over_represented_pvalue, method = "BH")
          
          GO.wall$under_represented_pvalue <-
            p.adjust(GO.wall$under_represented_pvalue, method = "BH")
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
        
        return(GO.wall)
      },
      error=function(cond) {
        
        sendSweetAlert(
          session = session,
          title = "Incorrect Gene type",
          text = "Ensure that the correct Gene symbol or genome is chosen",
          type = "error"
        )
        
        print(cond)
        
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
#' @return Returns a top 10 GO terms histogram
histGoTerms <- function(data, gof, session) {
  plot <- tryCatch({
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
        gsub("(.{18,}?)\\s", "\\1\n", levels(data$category))
    }

    plot <-
      ggplot(data, aes(reorder(category, log2(over_represented_pvalue)), -log2(over_represented_pvalue))) +
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
  },
  error = function(cond) {
    sendSweetAlert(
      session = session,
      title = "GO:Term Data Not Found",
      text = "Please run GO:TERM pipeline first",
      type = "error"
    )
    return()
  })

  return(plot)
}



#' Generate GO Term Info
#'
#' Calls GOTERM function and concatinates output
#'
#' @param goID ID of the GO term of interest
#' @param session Current R session
#'
#'
#' @return Returns the concatinated GO TERM information
goInfo <- function(goID, session) {
  out <- tryCatch({
    go.term <- GOTERM[[goID]]

    out <-  paste(GOID(go.term),
                  Term(go.term),
                  Definition(go.term),
                  Synonym(go.term),
                  sep = "\n")
  },
  error = function(cond) {
    sendSweetAlert(
      session = session,
      title = "GO:ID Not Found",
      text = "Please ensure that the GO:ID was appropriately typed",
      type = "error"
    )
    return()
  })
  return(out)
}
