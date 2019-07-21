#' Single Cell Cluster Tab UI
#'
#' @export
#' @return None
sc_clustUI <- function(id) {
  ns <- NS(id)
  tagList(
    # Sidebar panel for inputs ----
    sidebarPanel(
      actionButton(ns("elbowButton"), "Scale and Check Dimensionality"),

      conditionalPanel(
        condition =  "input.elbowButton > 0",
        ns = ns,

        tags$hr(),

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
            value = 0.5
          )
        ),

        conditionalPanel(
          condition = "input.clusterPackage==2",
          ns = ns,

          numericInput(
            ns("sc3minCells"),
            label = "Minimum Cells per Cluster",
            min = 20,
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
            choices = list("tSNE" = "tSNE",
                           "DDRTree (Pseudotime)" = "DDRTree")
          ),

          numericInput(
            ns("monoLowerDetection"),
            label = "Lower Detection Limit",
            min = 0,
            max = 1,
            value = 0.1
          ),


          conditionalPanel(
            condition = "input.monoRedMethod=='tSNE' ",
            ns = ns,

            numericInput(
              ns("monodimensionNo"),
              label = "Dimensions to be used",
              min = 3,
              value = 5
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

          h4("Visualize Data"),

          radioButtons(
            ns("clustplotType"),
            label = "Choose Plot Type",
            c(
              "Elbow plot" = "elbow",
              "PCA plot" = "pca",
              "tSNE Plot" = "tsne",
              "Dimensions Heatmap" = "heatmap"
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

          actionButton(ns("pcaButton"), "Generate Plot")

        )
      )
    ),

    # Main panel for displaying outputs ----
    mainPanel(

      plotOutput(ns("clusterPlot")),

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
#' @return Reactive value containing seurat object with scaled counts and reduced dimensions (PCA data)
sc_clust <- function(input, output, session, normData) {
  clust <- reactiveValues()

  observeEvent(input$elbowButton, {
    # if(!is.null(normData$normalizedData)){

    show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))
    clust$scaledData <- seuratElbow(normData$normalizedData)

    clust$clustPlot <- ElbowPlot(clust$scaledData)

    output$clusterPlot <- renderPlot({
      clust$clustPlot

    })

    hide_waiter()
  })


  observeEvent(input$clustButton, {
    if (!is.null(clust$scaledData)) {

      show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))
      if (input$clusterPackage == 1) {
        clust$finalData <-
          sueratClust(
            clust$scaledData,
            input$dimensionsInput,
            input$resolutionInput,
            as.numeric(input$clusterAlgCombo)
          )

      } else if (input$clusterPackage == 2) {
        print(1)
        clust$finalData <-
          sc3Cluster(
            clust$scaledData,
            input$estKCheck,
            input$sc3ClustNoInput,
            input$sc3minCells
          )

      } else{
        # Monocle normalization is suggested by the authors
        # however using Seurat Normalization seemed to work fine

        clust$finalData <-
          clusterMonocle(
            clust$scaledData,
            input$monoLowerDetection,
            input$monoRedMethod,
            input$monodimensionNo,
            input$estMonoClust,
            input$monoClustNo
          )

      }

      output$clustNoText <- renderText({
        sprintf("Number of estimated clusters is: %s",
                nlevels(clust$finalData@active.ident))
      })

    }

    hide_waiter()
  })



  observeEvent(input$pcaButton, {
    if (!is.null(clust$finalData)) {

      show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))

      if (is.null(clust$finalData@reductions$tsne)) {
        clust$finalData <- RunTSNE(clust$finalData)
      }

      if (input$clustplotType == "elbow") {
        clust$clustPlot <- ElbowPlot(clust$scaledData)
      } else if (input$clustplotType != "heatmap") {
        clust$clustPlot <-
          DimPlot(clust$finalData,
                  reduction = input$clustplotType,
                  pt.size = 1.2)
        clustPlotName <-
          paste("figures/", input$clustplotType, ".png", sep = "")
      } else{
        clust$clustPlot <-
          DimHeatmap(
            clust$finalData,
            dims = as.numeric(input$dimNoInput),
            cells = 1000,
            balanced = TRUE,
            fast = FALSE
          )
        clustPlotName <-
          paste("figures/",
                input$clustplotType,
                "DimNo",
                input$dimNoInput,
                ".png",
                sep = "")
      }

      output$clusterPlot <- renderPlot({
        clust$clustPlot

      })

      hide_waiter()
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
#' @param s_object Seurat object with normalized data
#'
#' @export
#' @return Seurat object with scaled counts and reduced dimensions (PCA data)
seuratElbow <- function(s_object) {
  scaled_object <- ScaleData(s_object)
  scaled_object <-
    RunPCA(scaled_object, features = VariableFeatures(object = scaled_object))

  return(scaled_object)
}

#' Seurat Clustering Pipeline
#'
#' @param s_object Seurat object with scaled and PCA data
#'
#' @export
#' @return Seurat object with clustering data
sueratClust <- function(s_object, dimNo, resQuant, algorithmNo) {
  s_object <- FindNeighbors(s_object, dims = 1:dimNo)
  s_object <-
    FindClusters(s_object, resolution = resQuant, algorithm = algorithmNo)

  return(s_object)
}

#' SC3 Clustering Pipeline
#'
#' @param s_object Seurat object with scaled and PCA data
#' @param estK Boolean - whether to estimate cluster number or not
#' @param clustNo If estK is false - Give desired cluster number
#' @param cellsPerC Minimum cells per cluster
#'
#' @export
#' @return Seurat object with SC3 clustering data
sc3Cluster <- function(s_object, estK, clustNo, cellsPerC) {
  # Convert from Seurat to sc3
  sce <- as.SingleCellExperiment(s_object)
  rowData(sce)$feature_symbol <- rownames(s_object)

  counts(sce) <- as.matrix(counts(sce))
  # normcounts(sce) <- as.matrix((s_object@assays$RNA@data))
  # logcounts(sce) <- as.matrix(logcounts(sce))
  logcounts(sce) <- as.matrix((s_object@assays$RNA@data))


  # Cluster similar cells
  qclust <-
    scran::quickCluster(sce, min.size = cellsPerC, assay.type = "logcounts")
  print(qclust)

  ## Data Normalization
  sce <-
    scran::computeSumFactors(
      sce,
      sizes = 20,
      clusters = qclust,
      positive = TRUE,
      assay.type = "logcounts"
    )
  # sce <- scater::normalize(sce)

  browseVignettes("SC3")
  if (estK) {
    # estimate No of clusters

    sce <- sc3_estimate_k(sce)
    clustNo <- sce@metadata$sc3$k_estimation
  }

  print(clustNo)
  # Perform unsupervised clustering of the cells and produce plots.
  sce <- sc3(object = sce,
             ks = clustNo)

  ### assign clusters from sc3 to s_object
  sce@metadata$sc3$consensus[[1]]$silhouette[, 1]
  clusters <- sce@metadata$sc3$consensus[[1]]$silhouette[, 1]
  names(clusters) <- colnames(s_object)
  s_object@active.ident <- as.factor(clusters)

  return(s_object)
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
#' @return Seurat object with SC3 clustering data
clusterMonocle <-
  function(s_object,
           lowerDetection,
           redMethod,
           dimensionNo,
           estimateClust,
           clustNo) {
    #1. Extract data, phenotype data, and feature data from the SeuratObject
    data <- as(as.matrix(s_object@assays$RNA@data), 'sparseMatrix')

    pd <- new('AnnotatedDataFrame', data = s_object@meta.data)

    fData <-
      data.frame(gene_short_name = row.names(data),
                 row.names = row.names(data))
    fd <- new('AnnotatedDataFrame', data = fData)


    # 2. Construct monocle cds
    my_cds <- newCellDataSet(
      data,
      phenoData = pd,
      featureData = fd,
      lowerDetectionLimit = lowerDetection,
      #* lowerDetection limit
      expressionFamily = uninormal()
    ) #* data type option (norm done with seurat)

    ## 3. normalisation and variance estimation steps (used in the differential expression analyses later on)
    # my_cds <- estimateSizeFactors(my_cds)
    # my_cds <- estimateDispersions(my_cds)

    # Not really a fix to DDRTree though..
    my_cds@reducedDimA <-
      t(s_object@reductions$pca@feature.loadings)


    ## 4. Dimension reduction
    my_cds <- reduceDimension(
      my_cds,
      max_components = 2,
      num_dim = dimensionNo,
      reduction_method = redMethod,
      scaling = TRUE,
      pseudo_expr = 0,
      norm_method = 'none',
      #Normalization from Seurat
      verbose = TRUE
    )

    if (redMethod == "tSNE") {
      if (estimateClust) {
        # 5A. Unsupervized Clustering
        my_cds <- clusterCells(my_cds)
      } else {
        # 5B. Unsupervised clustering requesting x-1 clusters
        my_cds <- clusterCells(my_cds, num_clusters = (clustNo + 1))
      }


      # 6. Store clusters
      s_object@active.ident <- pData(my_cds)$Cluster
    } else{
      # Get the "State" of each cell according to pseudotime
      my_cds <- orderCells(my_cds)

      # use the state as cluster in the seurat object
      s_object@active.ident <- pData(my_cds)$State

    }

    return(s_object)
  }
