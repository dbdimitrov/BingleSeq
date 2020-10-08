#' Differential Expression Tab UI
#'
#' @export
#' @return None
bulk_deDataUI <- function(id) {
  ns <- NS(id)

  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(

      h4("Run DE Pipeline"),

      selectInput(
        ns("selectPackage"),
        label = "Select DE package",
        choices = list(
          "DESeq2" = 1,
          "edgeR" = 2,
          "limma" = 3
        ),
        selected = 1
      ),

      conditionalPanel(
        condition = "input.selectPackage != 3",
        ns = ns,

        selectInput(
          ns("selectTestType"),
          label = "Select Test Type",
          choices = list(
            "Wald Test" = "Wald",
            "Likelihood Ratio Test" = "LRT"
          )
        )

      ),

      selectInput(
        ns("selectParamX"),
        label = "Select Fit Type",
        choices = list(
          "Parametric" = "parametric",
          "Local" = "local",
          "Mean" = "mean"
        )
      ),
      
      
      hr(),
      
      prettyCheckbox(ns("twoCondCheck"), label = tags$b("Use all samples for DE analysis"),
                    value = TRUE, bigger = TRUE, animation ="rotate"),
      
      conditionalPanel(condition = "!input.twoCondCheck",
                       ns = ns,
                       h5("Samples to be compared (2 condition contrast):"),
                       splitLayout(cellWidths = c("50%", "50%"), 
                                   uiOutput(ns("boxHolder_A")),
                                   uiOutput(ns("boxHolder_B"))
                                   )
      ),


      actionButton(ns("dePackageButton"), label = "Run DE Analysis"),

      conditionalPanel(condition = "input.dePackageButton > 0",
                       ns = ns,

                       hr(),
                       
                       h4("Filter DE results"),

                       checkboxInput(ns("exploreDE"),
                                     label = "",
                                     value = FALSE),

                       conditionalPanel(condition = "input.exploreDE",
                                        ns = ns,
                                        numericInput(
                                          ns("explorePvalue"),
                                          label = ("FDR threshold"),
                                          value = 0.05,
                                          min = 0.0001,
                                          max = 0.5
                                        ),
                                        textInput(ns("exploreGenes"),
                                                  "Please enter genes you wish to filer"),
                                        fluidRow(column(3, verbatimTextOutput(ns(
                                          "explorePvalue"
                                        )))),

                                        actionButton(ns("exploreButton"), label = "Filter Displayed Table")
                                        )
      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(# Output: DE Table ----
              htmlOutput(ns("helpDEInfo")),
              DT::dataTableOutput(ns("deTable")),
              conditionalPanel(condition = "input.dePackageButton > 0",
                               ns = ns,

                               downloadButton(ns("exploreDownload"), "Download Table")

              )
            )
  )
}


#' Differential Expression Tab Server
#'
#' @param fCounts The reactive value containing the filtered counts
#'
#' @export
#' @return Returns a reactive value with differential expression results as a list.
#' The first element of the list is the DE results
#' The second element is the normalized counts data
#' The third element is the metadata table
bulk_deData <- function(input, output, session, fCounts, unfCounts) {
  de <- reactiveValues()

  output$helpDEInfo <- renderUI({
    if(is.null(de$deTable)){
      HTML("<div style='border:2px solid blue; font-size: 14px; border-radius: 10px;'>

      <p style ='font-size: 15px; text-align: center'>
      <b>This tab enables DE analysis using <i>DESeq2,
      edgeR,</i> and <i>limma</i> pipelines. </b> </p> <br>

      <i> Please refer to the specific DE package manual
           for further information regarding the available
           test type, normalization methods, and fit type options. </i></div>")
    } else{
      HTML("<h4 style='padding-top: 8px'>Preview DE results</h4>
             <p style='padding-bot: 8px;'><i>
              The resulting parameters represent: <br>
              log2 expression fold-change (logFC);
              package/test specific test statistics;
              uncorrect p-value (Pvalue);
              multiple-testing FDR adjusted p-value (FDR).</i></p>")
    }
  })
  
  
  # Two Condition option ----
  sample_names <- reactive({
    #get names from metatable
    unfCounts$metaTable$sample
  })
  
  output$boxHolder_A <-
    renderUI({
      if (is.null(sample_names()))
        return()
      else
        checkboxGroupInput(session$ns("checkBox_A"), "Group A Samples", sample_names())
    })
  
  observeEvent(sample_names(), {
    updateCheckboxGroupInput(session, "checkBox_A", choices = sample_names())
  })
  
  observeEvent(input$checkBox_A, {
    de$treatment_A <- strsplit(input$checkBox_A, " ")
  })
  
  
  output$boxHolder_B <-
    renderUI({
      if (is.null(sample_names()))
        return()
      else
        checkboxGroupInput(session$ns("checkBox_B"), "Group B Samples", sample_names())
    })
  
  observeEvent(sample_names(), {
    updateCheckboxGroupInput(session, "checkBox_B", choices = sample_names())
  })
  
  observeEvent(input$checkBox_B, {
    de$treatment_B <- strsplit(input$checkBox_B, " ")
  })
  

  # Run DE ----
  observeEvent(input$dePackageButton, {
    # Save values to variables
    
    ## All or just two cond. -----
    if(!input$twoCondCheck){
        if(length(de$treatment_A) > 1 && length(de$treatment_B) > 1){
          df_sample_A <-
            data.frame(matrix(
              unlist(de$treatment_A),
              nrow = length(de$treatment_A),
              byrow = TRUE
            ))
          df_sample_A$batch <- "1"
          df_sample_A$treatment <- "A"
          
          df_sample_B <-
            data.frame(matrix(
              unlist(de$treatment_B),
              nrow = length(de$treatment_B),
              byrow = T
            ))
          df_sample_B$batch <- "1"
          df_sample_B$treatment <- "B"
          
          colnames(df_sample_A) <-
            c("sample", "batch", "treatment")
          colnames(df_sample_B) <-
            c("sample", "batch", "treatment")
          
          new_meta <- rbind(df_sample_A, df_sample_B)
          de$meta <- new_meta
          
          
          
          if(!is.null(fCounts$batchCorrected)){
            samples.used <- names(fCounts$batchCorrected)[(names(fCounts$batchCorrected) %in% new_meta$sample)]
            de$counts_batch <- fCounts$batchCorrected[, samples.used]
          } else{
            samples.used <- names(fCounts$filteredCounts)[(names(fCounts$filteredCounts) %in% new_meta$sample)]
            samples.used <- c(names(fCounts$filteredCounts)[1],samples.used)
            de$counts_filter <- fCounts$filteredCounts[, samples.used]
          }
        }
    } else {
      
      
      if(!is.null(fCounts$batchCorrected)){
        de$counts_batch <- fCounts$batchCorrected
      } else{
        de$counts_filter <- fCounts$filteredCounts
      }
      de$meta <- unfCounts$metaTable
      
    }
    
    
    
    de$selectedPackage <- as.numeric(input$selectPackage)


    show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))
      
    
    if(!is.null(fCounts$batchCorrected)){

      de$batched <- TRUE

      de$deTable <- dePipelineCaller(de$counts_batch,
                                     de$meta,
                                     input$selectTestType,
                                     input$selectParamX,
                                     as.numeric(input$selectPackage),
                                     session,
                                     TRUE
      )

    } else{

      de$batched <- FALSE
      de$deTable <- dePipelineCaller(de$counts_filter,
                                     de$meta,
                                     input$selectTestType,
                                     input$selectParamX,
                                     as.numeric(input$selectPackage),
                                     session,
                                     FALSE
      )

      if(!is.null(de$deTable))
        de$merged <- cbind(de$deTable[[1]],de$deTable[[2]])
    }
    

    observe({

      output$deTable <-
        DT::renderDataTable(
          de$deTable[[1]] %>%  rownames_to_column("gene_id") %>% 
            datatable(rownames = FALSE),
          options = list(pageLength = 10)
        )
    })

    hide_waiter()

    observeEvent(input$exploreButton,{
      de$deTable[[1]] <- de$deTable[[1]] %>%
        dplyr::filter(FDR < input$explorePvalue) %>%
        dplyr::filter(!grepl(as.character(input$exploreGenes), rownames(.)))

      output$deTable <-
        DT::renderDataTable(
          de$deTable[[1]] %>% rownames_to_column("gene_id") %>% 
            datatable(rownames = FALSE),
          options = list(pageLength = 10)
        )

    })

    output$exploreDownload <- downloadHandler(
      filename = function() {
        paste(format(Sys.time(), "%y-%m-%d_%H-%M"), "_deResults" , ".csv", sep = "")
      },
      content = function(file) {
        data <- de$merged

        write.csv(data, file)
      }
    )
  })

  # change UI choiceBoxes
  observeEvent(input$selectPackage,{
    if(as.numeric(input$selectPackage) == 1){
      updateSelectInput(session,
                        "selectTestType",
                        choices = list(
                          "Wald Test" = "Wald",
                          "Likelihood Ratio Test" = "LRT"
                        ))

      updateSelectInput(session,
                        "selectParamX",
                        label = "Select Fit Type",
                        choices = list(
                          "Parametric" = "parametric",
                          "Local" = "local",
                          "Mean" = "mean"
                        ))


    }else if(as.numeric(input$selectPackage) != 1){
      updateSelectInput(session,
                        "selectTestType",
                        choices = c("Exact Test" = "exact",
                                    "General Linear Model" = "GLM"))

      updateSelectInput(session,
                        "selectParamX",
                        label = "Select Normalization Method",
                        choices = list("Trimmed mean of M values(TMM)" = "TMM",
                                      "TMM with singleton pairing" = "TMMwsp",
                                      "Relative log expression" = "RLE",
                                      "Upper-quartile method" = "upperquartile",
                                      "No normalization" = "none"))

    }
  })

  return(de)

}


#' Differential Expression Pipeline Caller
#'
#' Function that calls the appropriate DGE pipeline
#' and provides exception handling
#'
#' @param readCounts The filtered CountTable
#' @param testType The type of test to be used
#' @param paramX Fit type for DESeq2, and normalization method for limma and edgeR
#' @param selectedPackage The package(pipeline) to be executed
#' @param session Current R session
#'
#'
#' @export
#' @return Returns the results of the appropriate DE pipeline
dePipelineCaller <- function(filteredCounts,
                             metaData,
                             testType,
                             paramX,
                             selectedPackage,
                             session,
                             useBatch) {
  out <- tryCatch(
    {
      if (selectedPackage == 1) {

        out <-
          deSequence(
            filteredCounts,
            metaData,
            testType,
            paramX,
            useBatch
          )

      } else if (selectedPackage == 2) {
        out <-
          deEdgeR(
            filteredCounts,
            metaData,
            testType,
            paramX,
            useBatch
          )

      } else if (selectedPackage == 3) {
        out <-
          deLimma(
            filteredCounts,
            metaData,
            paramX,
            useBatch
          )
      }
    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data Format Error",
        text = "Please ensure that the count and metadata tables were formatted and loaded correctly",
        type = "error"
      )
      return()
    }
  )
  return(out)
}



#' DESeq2 Differential Expression Pipeline
#'
#' @param readCounts The filtered CountTable
#' @param meta Metadata table
#' @param testType The type of test to be used
#' @param fitType The type of fitType to be used
#' @param useBatch If FALSE reasigns IDs to row.names
#'
#' @export
#' @return Returns DESeq2 DE pipeline results
deSequence <- function(readCounts, meta, testType, fitType, useBatch){

  if(!useBatch){
    output2 <- readCounts[, -1]
    rownames(output2) <- readCounts[, 1]
  } else{
    output2 <- readCounts
  }

  samplescondition <- meta$treatment

  if(!is.integer(output2[1,2])){ # if !int round to nearest int
    output2 <- apply(output2, 2, function(x) as.integer(round(x)))
  }

  coldata <-
    data.frame(row.names = colnames(output2), samplescondition)

  dds <-
    DESeqDataSetFromMatrix(
      countData = output2,
      colData = coldata,
      design =  ~ samplescondition
    )


  if(testType == "Wald"){
    dds <- DESeq(dds, test = "Wald", fitType = fitType)
  } else if(testType == "LRT"){
    dds = DESeq(dds, test = "LRT", reduced =  ~ 1, fitType = fitType)
  }

  res <-  results(dds)
  DEData <- data.frame(as.data.frame(res))

  DEData.matching <- DEData[, c(2, 4:6)]

  colnames(DEData.matching) <- c("logFC", as.character(testType), "Pvalue", "FDR")

  normCounts <- as.data.frame(counts(dds, normalized = TRUE))

  deList <- list(DEData.matching, normCounts, meta)

  return(deList)
}


#' EdgeR Differential Expression Pipeline
#'
#' @param readCounts The filtered CountTable
#' @param meta Metadata table
#' @param testType The type of test to be used
#' @param fitType The test method to be used
#' @param useBatch If FALSE reasigns IDs to row.names
#'
#' @export
#' @return EdgeR DE pipeline results
deEdgeR <- function(readCounts, meta, testType, normMethod, useBatch){
  ## 1. Load and Format Count dataframe (if !using batchCorr data)
  if(!useBatch){
    countTable  <- readCounts[, -1]
    rownames(countTable) <- readCounts[, 1]
  } else{
    countTable<- readCounts
  }

  ## 2. Assign Condition Groups
  samplescondition <- factor(meta$treatment)

  ## 3. Norm and Disp
  dge  <- DGEList(countTable, group = samplescondition)

  # Calculate normalization factors to scale the raw library sizes
  dge = calcNormFactors(dge, method = normMethod)

  # Estimate Common and Tagwise Dispersion
  dge = estimateCommonDisp(dge)
  dge = estimateTagwiseDisp(dge)

  ## 4. DE
  if(testType == "exact"){
    de = exactTest(dge)
  }else if(testType == "GLM"){
    fit <- glmQLFit(dge)
    de <- glmQLFTest(fit)
  }

  # 5. Extract DE results
  tt = topTags(de, n = nrow(dge), adjust.method = "fdr", sort.by	= "none")
  print(head(tt))
  if(ncol(tt$table)>4){
    res <- tt$table[,c(3,1,4,5)]
  } else{
    res <- tt$table[,c(2,1,4,3)]
  }

  ## 6.Extract normalized  CPMs
  nc = cpm(dge, normalized.lib.sizes = TRUE)
  normCounts <- as.data.frame(nc)

  ## 7. List with DE results and normalized counts
  deList <- list(res, normCounts, meta)

  return(deList)
}



#' Limma Differential Expression Pipeline
#'
#' @param readCounts The filtered CountTable
#' @param meta Metadata table
#' @param fitType The test method to be used
#' @param useBatch If FALSE reasigns IDs to row.names
#'
#' @export
#' @return Returns limma DE pipeline results
deLimma <- function(readCounts, meta, normMethod, useBatch) {

  if(!useBatch){
    countTable  <- readCounts[, -1]
    rownames(countTable) <- readCounts[, 1]
  } else{
    countTable<- readCounts
  }

  samplescondition <- factor(meta$treatment)
  d  <- DGEList(countTable, group = samplescondition)

  d <- calcNormFactors(d, method = normMethod)
  v <- voom(d, plot = FALSE)

  fit <- lmFit(v)
  ebayes.fit <- eBayes(fit)

  if(length(levels(as.factor(meta$treatment))) == 2){
  tab <-
    topTable(
      ebayes.fit,
      coef = 2,
      number = dim(ebayes.fit)[1],
      genelist = fit$genes$NAME,
      adjust = "fdr",
      sort.by = "none"
    )

  res <-data.frame(tab$logFC, tab$t, tab$P.Value, tab$adj.P.Val)
  colnames(res) <- c("logFC", "t", "Pvalue", "FDR")

  }else{
    tab <-
      topTable(
        ebayes.fit,
        coef = 1:length(levels(as.factor(meta$treatment))),
        number = dim(ebayes.fit)[1],
        genelist = fit$genes$NAME,
        adjust = "BH",
        sort.by = "none"
      )

    res <-data.frame(tab[,1], tab$F, tab$P.Value, tab$adj.P.Val)
    colnames(res) <- c("logFC", "F", "Pvalue", "FDR")
  }

  row.names(res) <- row.names(tab)

  deList <- list(res, v$E, meta)

  return(deList)
}
