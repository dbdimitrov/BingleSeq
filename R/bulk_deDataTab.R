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

      tags$hr(),

      numericInput(
        ns("conditionNo"),
        label = "Number of Conditions",
        value = 3,
        min = 2,
        max = 3
      ),

      hr(),
      fluidRow(column(3, verbatimTextOutput(
        ns("conditionNo")
      ))),

      numericInput(
        ns("replicateNo"),
        label = "Number of Replicates",
        value = 4,
        min = 1
      ),

      hr(),
      fluidRow(column(3, verbatimTextOutput(
        ns("replicateNo")
      ))),


      actionButton(ns("dePackageButton"), label = "Run DE Analysis"),

      conditionalPanel(condition = "input.dePackageButton > 0",
                       ns = ns,

                       checkboxInput(ns("exploreDE"), label = ("Explore DE results"),
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
#' @return Returns a reactive value with differential expression results
bulk_deData <- function(input, output, session, fCounts) {
  de <- reactiveValues()

  output$helpDEInfo <- renderUI({
    if(is.null(de$deTable)){
      HTML("<div style='border:2px solid blue; font-size: 14px; border-radius: 10px;'>

      <p style ='font-size: 15px; text-align: center'>
      <b>This tab enables DE analysis using <i>DESeq2,
      edgeR,</i> and <i>limma</i> pipelines. </b> </p> <br>

      Upon DE pipeline completion, a results table is displayed
      that contains the log2 expression fold-change (logFC),
      package specific test statistics, p-value,
      multiple-testing adjusted p-value (FDR). <br> <br>

      <i>Note: Diffrences among different samples are
      accounted for using the package-specific methods.
      Furthermore, only comparisons with up to 3
      conditions with the same number of replicates across
      each condition are implemented and the  appropriate
           number of conditions and replicates must be provided. </i></div>")
    } else{
      HTML("")
    }
  })

  # Run DE ----
  observeEvent(input$dePackageButton, {
    # Save values to variables
    de$selectedPackage <- as.numeric(input$selectPackage)
    de$conditionNo <- as.numeric(input$conditionNo)
    de$replicateNo <- as.numeric(input$replicateNo)

    if (de$conditionNo == 2) {
      de$offset <- 4
    } else if (de$conditionNo == 3) {
      de$offset <- 7
    } else if (de$conditionNo == 4) {
      de$offset <- 10
    }

    show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))

    de$deTable <- dePipelineCaller(fCounts$filteredCounts,
                                   as.numeric(input$conditionNo),
                                   as.numeric(input$replicateNo),
                                   as.numeric(input$selectPackage),
                                   session
                                   )
    observe({

      output$deTable <-
        DT::renderDataTable(
          de$deTable[, 1:(de$offset)] %>% datatable() %>%
            formatSignif(columns = c(1:de$offset), digits = 4),
          options = list(pageLength = 10)
        )
    })

    hide_waiter()

    observeEvent(input$exploreButton,{
      de$dispTable <-
        subset(de$deTable, FDR < input$explorePvalue)

      output$deTable <-
        DT::renderDataTable(
          de$dispTable[, 1:(de$offset)] %>% datatable() %>%
            formatSignif(columns = c(1:de$offset), digits = 4),
          options = list(pageLength = 10)
        )

    })

    output$exploreDownload <- downloadHandler(
      filename = function() {
        paste(format(Sys.time(), "%y-%m-%d_%H-%M"), "_deResults" , ".csv", sep = "")
      },
      content = function(file) {
        data <- de$deTable

        write.csv(data, file)
      }
    )
  })

  return(de)

}


#' Differential Expression Pipeline Caller
#'
#' Function that calls the appropriate DGE pipeline
#' and provides exception handling
#'
#' @param readCounts The filtered CountTable
#' @param conditionNo The number of Conditions
#' @param replicateNo The number of Replicates
#' @param selectedPackage The package(pipeline) to be executed
#' @param session Current R session
#'
#'
#' @export
#' @return Returns the results of the appropriate DE pipeline
dePipelineCaller <- function(filteredCounts, conditionNo, replicateNo, selectedPackage, session) {
  out <- tryCatch(
    {
      if (selectedPackage == 1) {
        out <-
          deSequence(
            filteredCounts,
            conditionNo,
            replicateNo
          )

      } else if (selectedPackage == 2) {
        out <-
          deEdgeR(
            filteredCounts,
            conditionNo,
            replicateNo
          )

      } else if (selectedPackage == 3) {
        out <-
          deLimma(
            filteredCounts,
            conditionNo,
            replicateNo
          )
      }
    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Data Format Error",
        text = "Please ensure that the data was formatted and loaded correctly",
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
#' @param conditionNo The number of Conditions
#' @param replicateNo The number of Replicates
#'
#' @export
#' @return Returns DESeq2 DE pipeline results
deSequence <- function(readCounts, conditionNo, replicateNo) {

  output2 <- readCounts[, -1]
  rownames(output2) <- readCounts[, 1]

  if (conditionNo == 2) {
    samplescondition <- c(rep("C1", replicateNo), rep("C2", replicateNo))

  } else if (conditionNo == 3) {
    samplescondition <-
      c(rep("C1", replicateNo),
        rep("C2", replicateNo),
        rep("C3", replicateNo))
  }


  coldata <-
    data.frame(row.names = colnames(output2), samplescondition)


  dds <-
    DESeqDataSetFromMatrix(
      countData = output2,
      colData = coldata,
      design =  ~ samplescondition
    )


  if (conditionNo == 2) {
    dds <- DESeq(dds)
    res <-
      results(dds, contrast = c("samplescondition", "C2", "C1"))

    DEData <- data.frame(as.data.frame(res))


    DEData.matching <- DEData[, c(2, 4:6)]

    colnames(DEData.matching) <- c("logFC", "stat", "Pvalue", "FDR")

  } else{
    dds = DESeq(dds, test = "LRT", reduced =  ~ 1)

    res <- results(dds)

    AvB <-
      results(dds,
              contrast = c("samplescondition", "C2", "C1"),
              test = "Wald")

    BvC <-
      results(dds,
              contrast = c("samplescondition", "C3", "C2"),
              test = "Wald")

    res$log2FC_AvB <- AvB$log2FoldChange
    res$log2FC_BvC <- BvC$log2FoldChange


    DEData.matching <- res[, c(7, 8, 2:6)]
    colnames(DEData.matching) <-
      c("logFC_C2vC1",
        "logFC_C3vC2",
        "logFC_C3vC1",
        "lfcSE",
        "stat",
        "Pvalue",
        "FDR")

    DEData.matching <- data.frame(as.data.frame(DEData.matching))

  }

  resdata <- as.data.frame(counts(dds, normalized = TRUE))

  masterFileDE <-
    merge(DEData.matching, resdata, by = "row.names", all.x = TRUE)

  format.data.frame(masterFileDE, digits = 4)

  countTable <- masterFileDE[, -1]
  rownames(countTable) <- masterFileDE[, 1]

  return(countTable)
}

#' EdgeR Differential Expression Pipeline
#'
#' @param readCounts The filtered CountTable
#' @param conditionNo The number of Conditions
#' @param replicateNo The number of Replicates
#'
#' @export
#' @return EdgeR DE pipeline results
deEdgeR <- function(readCounts, conditionNo, replicateNo) {
  ## 1. Load and Format Count dataframe
  countTable  <- readCounts[, -1]
  rownames(countTable) <- readCounts[, 1]


  # 2. Create condition groups
  if (conditionNo == 2) {
    samplescondition <- c(rep("C1", replicateNo), rep("C2", replicateNo))

  } else if (conditionNo == 3) {
    samplescondition = factor(c(
      rep("C1", replicateNo),
      rep("C2", replicateNo),
      rep("C3", replicateNo)
    ))
    my.contrasts <-
      makeContrasts(
        BvsA = C2 - C1,
        CvsB = C3 - C2,
        CvsA = C3 - C1,
        levels = samplescondition
      )
    design <- model.matrix( ~ 0 + samplescondition) # without intercept
    colnames(design) <- levels(samplescondition)
  }


  ## 3. Norm and DE
  dge  <- DGEList(countTable, group = samplescondition)

  # Calculate normalization factors to scale the raw library sizes
  dge = calcNormFactors(dge)

  # Estimate Common and Tagwise Dispersion
  dge = estimateCommonDisp(dge)
  dge = estimateTagwiseDisp(dge)

  if (conditionNo == 2) {
    de = exactTest(dge, pair = c("C1", "C2"))
    tt = topTags(de, n = nrow(dge))

  } else {
    fit <- glmFit(dge, design)
    lrt <- glmLRT(fit, contrast = my.contrasts)

    tt = topTags(lrt, n = nrow(dge))
    head(tt$table)
  }

  ## 4.Extract normalized  CPMs
  nc = cpm(dge, normalized.lib.sizes = TRUE)
  resdata <- as.data.frame(nc)


  ## 5. Format and create MF
  masterFileDE <- merge(tt$table, resdata, by = "row.names", all.x = TRUE)

  countTable <- masterFileDE[, -1]
  rownames(countTable) <- masterFileDE[, 1]

  return(countTable)
}



#' Limma Differential Expression Pipeline
#'
#' @param readCounts The filtered CountTable
#' @param conditionNo The number of Conditions
#' @param replicateNo The number of Replicates
#'
#' @export
#' @return Returns limma DE pipeline results
deLimma <- function(readCounts, conditionNo, replicateNo) {

  countTable  <- readCounts[, -1]
  rownames(countTable) <- readCounts[, 1]

  # Create Name vectors for conditions and contrasts
  if (conditionNo == 2) {
    samplescondition <- c(rep("C1", replicateNo), rep("C2", replicateNo))
    design = model.matrix(object = ~ samplescondition)
    coef = 2

  } else if (conditionNo == 3) {
    samplescondition = factor(c(
      rep("C1", replicateNo),
      rep("C2", replicateNo),
      rep("C3", replicateNo)
    ))
    my.contrasts <-
      makeContrasts(
        BvsA = C2 - C1,
        CvsB = C3 - C2,
        CvsA = C3 - C1,
        levels = samplescondition
      )
    design <- model.matrix( ~ 0 + samplescondition)
    colnames(design) <- levels(samplescondition)
    coef = c(1:3)
  }

  d  <- DGEList(countTable, group = samplescondition)

  # TMM normalize
  d <- calcNormFactors(d)

  # Transform count data to logCPM
  v <- voom(d, design = design, plot = FALSE)

  fit <- lmFit(v, design)

  if (conditionNo > 2) {
    fit2 <- contrasts.fit(fit, my.contrasts)
    ebayes.fit <- eBayes(fit2)
  } else{
    ebayes.fit <- eBayes(fit)
  }

  tab <-
    topTable(
      ebayes.fit,
      coef = coef ,
      number = dim(ebayes.fit)[1],
      genelist = fit$genes$NAME,
      adjust = "fdr",
      sort.by = "none"
    )

  head(tab, n = 100)

  if (conditionNo == 2) {
    DEData <- tab[, c(1, 3:5)]

    colnames(DEData)[4] <- "FDR"
  } else if (conditionNo == 3) {
    DEData <- tab[, c(1:7)]
    colnames(DEData)[7] <- "FDR"
  }

  IDs <- rownames(tab)
  masterFileDE <- cbind(IDs, DEData, v$E)
  rownames(masterFileDE) <- NULL

  countTable <- masterFileDE[, -1]
  rownames(countTable) <- readCounts[, 1]

  return(countTable)
}
