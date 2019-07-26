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
            ns("sc3minDropout"),
            label = "Filter genes by minimum percent of dropouts",
            min = 0,
            max = 50,
            value = 0
          ),

          numericInput(
            ns("sc3maxDropout"),
            label = "Filter genes by maximum percent of dropouts",
            min = 1,
            max = 100,
            value = 100
          ),

          numericInput(
            ns("sc3nStart"),
            label = "Random sets used in clustering (nStart)",
            min = 50,
            value = 1000
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

          h4("Visualize Data"),

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
#' @return Returns a Reactive value containing seurat object with scaled counts and reduced dimensions (PCA data)
sc_clust <- function(input, output, session, normData) {
  clust <- reactiveValues()

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
        # Monocle normalization is suggested by the authors
        # however using Seurat Normalization seemed to work fine

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

      output$clustNoText <- renderText({
        sprintf("Number of estimated clusters is: %s",
                nlevels(clust$finalData@active.ident))
      })

    }

    hide_waiter()
  })



  observeEvent(input$pcaButton, {
    if (!is.null(clust$finalData)) {

      # show_waiter(tagList(spin_folding_cube(), h2("Loading ...")))

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

      # hide_waiter()
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
  scaled_object <- ScaleData(s_object)
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
                      pt.size = 1.4)


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
      # # delete
      # s_object <- pbmc
      # nStart = 1000
      # minDrop = 0
      # maxDrop = 100
      # head(counts(sce))
      # head(normcounts(sce))


      # Convert sparse matrix counts Seurat object to dense matrix in SC3 object
      sce <- as.SingleCellExperiment(s_object)
      rowData(sce)$feature_symbol <- rownames(s_object)

      counts(sce) <- as.matrix(counts(sce)) # using non-normalized counts
      logcounts(sce) <- as.matrix((s_object@assays$RNA@data)) # normalized counts (used in clustering)


      # Delete commented lines (only used in testing)
      # qclust <-
      #   scran::quickCluster(sce, min.size = 10, assay.type = "logcounts")
      # print(qclust)


      # sce <-
      #   scran::computeSumFactors(
      #     sce,
      #     sizes = 20,
      #     clusters = qclust,
      #     positive = TRUE,
      #     min.mean = 2, # NumericInput required (0.5 as default) + tryCatch
      #     assay.type = "logcounts"
      #   )

      # Data Normalization (done with Seurat)
      # sce <- scater::normalize(sce) #probvam s i bez

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
                 ks = clustNo)


      ### assign clusters from sc3 to s_object
      sce@metadata$sc3$consensus[[1]]$silhouette[, 1]
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
#' @return Returns a list containing a Seurat object with Monocle-produced clustering data and tSNE plot
clusterMonocle <-
  function(s_object,
           lowerDetection,
           redMethod,
           clustMethod,
           dimensionNo,
           estimateClust,
           clustNo,
           session) {

      # delete
    # lowerDetection = 0.1
    # dimensionNo = 10
    # redMethod = "tSNE"
    # clustMethod = "densityPeak"
    # clustNo = 4

    out <- tryCatch(
      {
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

        # Save PCA to this object
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
            my_cds <- clusterCells(my_cds, method = clustMethod)
          } else {
            # 5B. Unsupervised clustering requesting x-1 clusters
            my_cds <- clusterCells(my_cds, num_clusters = (clustNo + 1), method = clustMethod)
          }

          # 6. Store clusters

          print(pData(my_cds)$Cluster)

          clusters <- pData(my_cds)$Cluster
          names(clusters) <- colnames(s_object)
          s_object@active.ident <- as.factor(clusters)


          tsne <- plot_cell_clusters(my_cds) # Did not work with DDRTree

          out <- list(s_object, tsne)

        } else{
          # Get the "State" of each cell according to pseudotime
          my_cds <- orderCells(my_cds)

          # use the state as cluster in the seurat object
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
