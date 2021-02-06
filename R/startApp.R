#' Launch the application
#'
#' #' Launches the app
#'
#' @import shiny
#' @import shinyjs
#' @import shinyFiles
#' @import shinyWidgets
#' @import waiter
#' @import DT
#' @import grDevices
#' @import ggplot2
#' @import VennDiagram
#' @import factoextra
#' @import fastcluster
#' @import gridExtra
#' @import ggrepel
#' @import GenomeInfoDbData
#' @import geneLenDataBase
#' @import sva
#' @import Harman
#' @import DESeq2
#' @import edgeR
#' @import limma
#' @import MAST
#' @import monocle
#' @import Seurat
#' @import scran
#' @import SC3
#' 
#' @import goseq
#' @import GO.db
#' 
#' 
#' @import pheatmap
#' @import tidyverse
#' @import progeny
#' @import dorothea
#' @import plotly
#' 
#' @include server.R
#' @include ui.R
#'
#' @export
#' @return None
startBingleSeq <- function() {
  if (interactive()) {
    # Load packages
    require(shiny)
    require(shinyjs)
    require(shinyFiles)
    require(shinyWidgets)

    require(waiter)
    require(dplyr)
    require(DT)
    require(grDevices)

    require(ggplot2)
    require(VennDiagram)
    require(factoextra)
    require(fastcluster)
    require(reshape2)
    require(gridExtra)
    require(ggrepel)

    require(Harman)
    require(sva)

    require(DESeq2)
    require(edgeR)
    require(limma)
    require("MAST")

    require(monocle)
    require(Seurat)
    require('scran')
    require('SC3')

    require(goseq)
    require(GO.db)

    require(pheatmap)
    require(tidyverse)
    require(progeny)
    require(dorothea)
    require(plotly)
    require(viper)
    
    options(shiny.maxRequestSize=500*1024^2)

    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)
  }
}
