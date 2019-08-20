#' Single Cell Cluster Tab UI
#'
#' @export
#' @return None
sc_clustUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(

      h4("Perform Unsupervised clustering"),

      actionButton(ns("elbowButton"), "Generate Clustering Prerequisites"),

      conditionalPanel(
        condition =  "input.elbowButton > 0",
        ns = ns,

        tags$hr(),

        h4("Run Clustering"),

        selectInput(
          ns("clusterPackage"),
          label = "Select Clustering Package",
          choices = list(
            "Seurat" = 1,
            "SC3" = 2,
            "Monocle" = 3
          )
        ),

        conditionalPanel(
          condition = "input.clusterPackage==1",
          ns = ns,

          selectInput(
            ns("clusterAlgCombo"),
            label = "Select Clustering method",
            choices =
              list(
                "Louvain algorithm" = 1,
                "Louvain algorithm with multilevel refinement" = 2,
                "SLM algorithm" = 3
              )
          ),

          numericInput(
            ns("dimensionsInput"),
            label = "Dimensions to be used",
            min = 5,
            max = 50,
            value = 10
          ),

          numericInput(
            ns("resolutionInput"),
            label = "Resolution parameter Value",
            min = 0,
            max = 3,
            value = 0.6
          )
        ),

        conditionalPanel(
          condition = "input.clusterPackage==2",
          ns = ns,

          numericInput(
            ns("sc3minDropout"),
            label = "Filter genes below minimum percent of dropouts",
            min = 0,
            max = 50,
            value = 0
          ),

          numericInput(
            ns("sc3maxDropout"),
            label = "Filter genes above maximum percent of dropouts",
            min = 1,
            max = 100,
            value = 100
          ),

          numericInput(
            ns("sc3nStart"),
            label = "Random sets used in clustering (nStart)",
            min = 50,
            value = 50
          ),

          checkboxInput(ns("estKCheck"), label = "Estimate Cluster Number", TRUE),

          conditionalPanel(
            condition = "!input.estKCheck",
            ns = ns,

            numericInput(
              ns("sc3ClustNoInput"),
              label = "Number of Clusters",
              min = 3,
              value = 6
            )
          )
        ),

        conditionalPanel(
          condition = "input.clusterPackage==3",
          ns = ns,

          selectInput(
            ns("monoRedMethod"),
            label = "Select Reduction Method",
            choices = list("tSNE" = "tSNE"
                           # ,"DDRTree (Pseudotime)" = "DDRTree"
                           )
          ),

          conditionalPanel(
            condition = "input.monoRedMethod=='tSNE' ",
            ns = ns,

            selectInput(
              ns("monoClustMethod"),
              label = "Select Clustering Method",
              choices = list("Density Peak" = "densityPeak",
                             "Louvian algorithm" = "louvain")
            ),

            numericInput(
              ns("monodimensionNo"),
              label = "Dimensions to be used",
              min = 3,
              value = 10
            ),

          numericInput(
            ns("monoLowerDetection"),
            label = "Lower Detection Limit",
            min = 0,
            max = 1,
            value = 0
          ),


            checkboxInput(ns("estMonoClust"),
                          label = "Estimate Cluster Number", TRUE)
          ),


          conditionalPanel(
            condition = "!input.estMonoClust",
            ns = ns,

            numericInput(
              ns("monoClustNo"),
              label = "Number of Clusters",
              min = 2,
              max = 50,
              value = 5
            )
          )
        ),


        actionButton(ns("clustButton"), "Cluster Data"),

        conditionalPanel(
          condition =  "input.clustButton > 0",
          ns = ns,
          tags$hr(),

          h4("Visualize Results"),

          radioButtons(
            ns("clustplotType"),
            label = "Choose Plot Type",
            c(
              "Elbow plot" = "elbow",
              "Dimensions Heatmap" = "heatmap",
              "PCA plot" = "pca",
              "tSNE Plot" = "tsne"
            )
          ),

          conditionalPanel(
            condition = "input.clustplotType=='heatmap'",
            ns = ns,
            numericInput(
              ns("dimNoInput"),
              label = "Dimension Number",
              min = 1,
              max = 20,
              value = 1
            )
          ),

          actionButton(ns("pcaButton"), "Generate Plot"),

          tags$hr(),

          actionButton(ns("saveObjectButton"), "Save Object")
        )
      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      htmlOutput(ns("helpClustInfo")),


      plotOutput(ns("clusterPlot"), width = "800px", height = "500px"),

      verbatimTextOutput(ns("clustNoText"), placeholder = F),

      downloadButton(ns("downloadClustPlot"), "Download Plot")
    )

  )
}


#' Single Cell Cluster Tab Server
#'
#' @param normData Reactive value containing seurat object with normalized data
#'
#' @export
#' @return Returns a Reactive value containing seurat object
#'  with scaled counts and reduced dimensions (PCA data)
sc_clust <- function(input, output, session, normData) {
  clust <- reactiveValues()

  output$helpClustInfo <- renderUI({
    if(input$elbowButton == 0){
    HTML("<div style='border:2px solid blue; padding-top: 8px;
    padding-bottom: 8px; font-size: 14px; border-radius: 10px;'>
     <p style='text-align: center'> <b>This tab enables unsupervised clustering with
     <i> Seurat, monocle, and SC3</i>. </b> </p> <br>

    To begin, press the 'Generate Clustering Prerequisites' button.
    What this will do is to: <br>
    Scale the data, linear dimensional reduction with PCA,
    and return an Elbow plot that can be used in PC Selection. <br> <br>

    Once the prerequisites are generated,
    the unsupervised clustering methods will become available. <br>
    Subsequent to clustering, visualization options themselves become available.
        </div> ")
    } else {
    if(input$clusterPackage == 1){
      HTML("<div style='border:2px solid blue; font-size: 14px;
       padding-top: 8px; padding-bot: 8px; border-radius: 10px;'>
       <p style='text-align: center'><b>Unsupervised Clustering with <i>Seurat</i>:</b> </p> <br>

       First, choose clustering algorithm of preference. <br>
       Then, provide the true dimensionality of the dataset via the 'Dimension to be used' parameter. <br>
       The true dimensionality can be estimated by looking at the 'elbow' of the elbow plot. <br> <br>

       Finally, use the 'Resolution' parameter to set the
       clustering ‘granularity’ and control the number of clusters.
       <br>Note that the optimal 'Resolution' for datasets with ~3000 cells is 0.6-1.2
       and it is typically higher for larger datasets. <br> </div>")
    } else if(input$clusterPackage == 2){
      HTML("<div style='border:2px solid blue; font-size: 14px;
      padding-top: 8px; padding-bot: 8px; border-radius: 10px;'>
      <p style='text-align: center'><b>Unsupervised Clustering with <i>SC3</i>:</b> </p> <br>

      SC3’s inbuilt filtering options enable further reducion of noise by filtering out <br>
      genes below and above certain dropout (zero value) percentage thresholds. <br> <br>

      'nStart' parameter enables control over the number of
      random datasets used in clustering and hence computational time. <br>
      By default, this parameter is set to 1000 when working with
      less than 2000 cells and to 50 when working with more than 2000 cells. <br>
      SC3 is magnitudes slower than the other approaches
      and appropriately setting 'nStart' is essential. <br> <br>

      Note that SC3 enables the number of clusters to be specified or estimated. <br>
      However, it should be noted that estimating
           cluster number with SC3 often results in overestimations. </div")

    } else{
      HTML(
        "<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bot: 8px; border-radius: 10px;'>
              <p style='text-align: center'><b>Unsupervised
              Clustering with <i>monocle</i>:</b> </p> <br>

        First, Select clustering algorithm of preference
        (use Louvian when working with large datasets). <br>
        Then use 'Dimension to be used' parameter, to provide the true dimensionality,
        determined using the Elbow plot <br>
        If required, further filter noise according
        to the 'Lower Detection Limit' parameter <br> <br>

        Finally, choose whether to estimate with
        monocle or specify a desired cluster number. </div>"
      )
    }
    }
  })


  observeEvent(input$elbowButton, {

    show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))

    clust$scaledData <- seuratElbow(normData$normalizedData)
    clust$clustPlot <- clust$scaledData[[2]]

    output$clusterPlot <- renderPlot({
      clust$clustPlot
    })

    hide_waiter()
  })


  observeEvent(input$clustButton, {
    if (!is.null(clust$scaledData)) {

      show_waiter(tagList(spin_folding_cube(), h2("Loading...")))


      if (input$clusterPackage == 1) {
        clust$results <-
          sueratClust(
            clust$scaledData[[1]],
            input$dimensionsInput,
            input$resolutionInput,
            as.numeric(input$clusterAlgCombo),
            session
          )

        clust$finalData <- clust$results[[1]]

      } else if (input$clusterPackage == 2) {

        clust$results <-
          sc3Cluster(
            clust$scaledData[[1]],
            input$sc3minDropout,
            input$sc3maxDropout,
            input$estKCheck,
            input$sc3ClustNoInput,
            input$sc3nStart,
            session
          )

        clust$finalData <- clust$results[[1]]

      } else{
        clust$results <-
          clusterMonocle(
            clust$scaledData[[1]],
            input$monoLowerDetection,
            input$monoRedMethod,
            input$monoClustMethod,
            input$monodimensionNo,
            input$estMonoClust,
            input$monoClustNo,
            session
          )

        clust$finalData <- clust$results[[1]]
      }

      clust$clustPlot <- clust$results[[2]]

      output$clustNoText <- renderText({
        sprintf("Number of estimated clusters is: %s",
                nlevels(clust$finalData@active.ident))
      })

    }

    hide_waiter()
  })


  observeEvent(input$pcaButton, {
    if (!is.null(clust$finalData)) {
      if (input$clustplotType == "elbow") {
        clust$clustPlot <- clust$scaledData[[2]]

      } else if (input$clustplotType == "pca") {
        clust$clustPlot <- DimPlot(clust$finalData,
                                   reduction = input$clustplotType,
                                   pt.size = 1.4)

      } else if(input$clustplotType == "tsne"){

        clust$clustPlot <- clust$results[[2]]

      } else{
        clust$clustPlot <-
          DimHeatmap(
            clust$finalData,
            dims = as.numeric(input$dimNoInput),
            cells = 1000,
            balanced = TRUE,
            fast = FALSE
          )

        hm.palette <-
          colorRampPalette(c("red", "white", "blue")) # Set the colour range

        clust$clustPlot <- clust$clustPlot +  scale_fill_gradientn(colours = hm.palette(100))
      }

      output$clusterPlot <- renderPlot({
        clust$clustPlot
      })
    }
  })

  # Download current plot
  output$downloadClustPlot <- downloadHandler(
    filename = function() {
      paste(as.character(input$clustplotType), device = ".png", sep = "")
    },
    content = function(file) {
      device <- function(..., width, height) {
        grDevices::png(
          ...,
          width = width,
          height = height,
          units = "px",
          pointsize = 14
        )
      }
      ggsave(
        file,
        plot = clust$clustPlot,
        device = device,
        width = 1280,
        height = 720,
        limitsize = FALSE
      )
    }
  )
  return(clust)
}


#' Single Cell Scale and Dimension Reduction function
#'
#' @param scaled_object Seurat object with scaled data
#'
#' @export
#' @return Returns a Seurat object with scaled counts and reduced dimensions (PCA data)
seuratElbow <- function(s_object) {

  all.genes <- rownames(s_object)
  scaled_object <- ScaleData(s_object, features = all.genes)
  scaled_object <-
    RunPCA(scaled_object, features = VariableFeatures(object = scaled_object))

  elbow <- ElbowPlot(scaled_object)

  out <- list(scaled_object, elbow)

  return(out)
}

#' Seurat Clustering Pipeline
#'
#' @param s_object Seurat object with scaled and PCA data
#' @param dimNo Number of dimensions to be used when clustering
#' @param resQuant Resolution parameter used to control the number of clusters
#' @param algorithmNo The clustering algorithm to be used
#' @param session Current R session
#'
#' @export
#' @return Returns a list containing a Seurat object with clustering data and a tSNE plot
sueratClust <- function(s_object, dimNo, resQuant, algorithmNo, session) {

  out <- tryCatch(
    {
      s_object <- FindNeighbors(s_object, dims = 1:dimNo)
      s_object <- FindClusters(s_object,
                               resolution = resQuant,
                               algorithm = algorithmNo)
      s_object <- RunTSNE(s_object, dims = 1:dimNo)


      tsne <- DimPlot(s_object,
                      reduction = "tsne",
                      pt.size = 1.6)

      out <- list(s_object, tsne)
    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Clustering Error Encountered",
        text = "Consider applying more stringent QC",
        type = "error"
      )
      return()
    }
  )

  return(out)
}

#' SC3 Clustering Pipeline
#'
#' @param s_object Seurat object with scaled and PCA data
#' @param minDrop Genes with percent of dropouts smaller than minDrop are filtered out before clustering.
#' @param maxDrop Genes with percent of dropouts larger than maxDrop are filtered out before clustering.
#' @param estK Boolean - whether to estimate cluster number or not
#' @param clustNo If estK is false - Give desired cluster number
#' @param nStart Minimum cells per cluster
#' @param Current R session
#'
#' @export
#' @return Returns a list containing a Seurat object with SC3-produced clustering data and tSNE plot
sc3Cluster <- function(s_object, minDrop, maxDrop, estK, clustNo, nStart, session) {

  out <- tryCatch(
    {

      # Convert sparse matrix counts Seurat object to dense matrix in SC3 object
      sce <- as.SingleCellExperiment(s_object)
      rowData(sce)$feature_symbol <- rownames(s_object)

      counts(sce) <- as.matrix(counts(sce))
      logcounts(sce) <- as.matrix((s_object@assays$RNA@data))


      if (estK) {
        sce <- sc3_estimate_k(sce) # estimate clustNo
        clustNo <- sce@metadata$sc3$k_estimation
      }

      # Perform unsupervised clustering of the cells and produce plots.
      sce <- sc3(object = sce,
                 gene_filter = TRUE,
                 pct_dropout_min = minDrop,
                 pct_dropout_max = maxDrop,
                 kmeans_nstart = nStart,
                 svm_max = 100000, # prevents reshuffling when working with large datasets
                 ks = clustNo)


      ### assign clusters from sc3 to s_object
      print(sce@metadata$sc3$consensus[[1]]$silhouette[, 1])
      clusters <- sce@metadata$sc3$consensus[[1]]$silhouette[, 1]
      names(clusters) <- colnames(s_object)
      s_object@active.ident <- as.factor(clusters)

      # tSNE plot
      k_estimated_field <- paste("sc3", clustNo, "clusters", sep = '_')

      tsne <- scater::plotTSNE(sce , colour_by = k_estimated_field) +
        theme_classic() +
        guides(fill=guide_legend("Clusters")) +
        theme(legend.text=element_text(size=12))

      out <- list(s_object, tsne)

    },
    error=function(cond) {
      sendSweetAlert(
        session = session,
        title = "Clustering Error Encountered",
        text = "Consider using another clustering method or applying more stringent QC",
        type = "error"
      )
      return()
    }
  )

  return(out)
}

#' Monocle Clustering Pipeline
#'
#' @param s_object Seurat object with scaled and PCA data
#' @param lowerDetection The minimum expression level that consistitutes true expression
#' @param redMethod Dimension reduction method to be used
#' @param dimensionNo Number of Dimensions to be used in clustering (tSNE)
#' @param estimateClust Boolean - Whether to estimate cluster No or not
#' @param clustNo If estimateClust is False - provide number of desired clusters
#'
#' @export
#' @return Returns a list containing a Seurat object with
#' Monocle-produced clustering data and tSNE plot
clusterMonocle <-
  function(s_object,
           lowerDetection,
           redMethod,
           clustMethod,
           dimensionNo,
           estimateClust,
           clustNo,
           session) {

    out <- tryCatch(
      {
        data <- as(as.matrix(s_object@assays$RNA@data), 'sparseMatrix')

        pd <- new('AnnotatedDataFrame', data = s_object@meta.data)

        fData <-
          data.frame(gene_short_name = row.names(data),
                     row.names = row.names(data))
        fd <- new('AnnotatedDataFrame', data = fData)

        my_cds <- newCellDataSet(
          data,
          phenoData = pd,
          featureData = fd,
          lowerDetectionLimit = lowerDetection,
          expressionFamily = uninormal()
        )

        my_cds@reducedDimA <-
          t(s_object@reductions$pca@feature.loadings)

        ## Dimension reduction
        my_cds <- reduceDimension(
          my_cds,
          max_components = 2,
          num_dim = dimensionNo,
          reduction_method = redMethod,
          scaling = TRUE,
          pseudo_expr = 0,
          norm_method = 'none',
          verbose = TRUE
        )

        if (redMethod == "tSNE") {
          if (estimateClust) {
            # Unsupervized Clustering
            my_cds <- clusterCells(my_cds, method = clustMethod)
          } else {
            # Unsupervised clustering requesting x-1 clusters
            my_cds <- clusterCells(my_cds, num_clusters = (clustNo + 1), method = clustMethod)
          }

          clusters <- pData(my_cds)$Cluster
          names(clusters) <- colnames(s_object)
          s_object@active.ident <- as.factor(clusters)

          tsne <- plot_cell_clusters(my_cds, cell_size = 1.6)

          out <- list(s_object, tsne)

        } else{
          my_cds <- orderCells(my_cds)
          s_object@active.ident <- pData(my_cds)$State

        }
      },
      error=function(cond) {
        sendSweetAlert(
          session = session,
          title = "Clustering Error Encountered",
          text = "Consider using another clustering method/package or applying more stringent QC",
          type = "error"
        )

        message(cond)
        return()
      }
    )
    return(out)
  }
