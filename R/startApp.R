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
#' @import reshape2
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
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import org.Dm.eg.db
#' @import org.Dr.eg.db
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
    library(shiny)
    library(shinyjs)
    library(shinyFiles)
    library(shinyWidgets)

    library(waiter)
    library(dplyr)
    library(DT)
    library(grDevices)

    library(ggplot2)
    library(VennDiagram)
    library(factoextra)
    library(fastcluster)
    library(reshape2)
    library(gridExtra)
    library(ggrepel)

    library(Harman)
    library(sva)

    library(DESeq2)
    library(edgeR)
    library(limma)
    library("MAST")

    library(monocle)
    library(Seurat)
    library('scran')
    library('SC3')
    
    library(AnnotationDbi)
    library(goseq)
    library(GO.db)

    library(pheatmap)
    library(tidyverse)
    library(progeny)
    library(dorothea)
    library(plotly)
    
    options(shiny.maxRequestSize=500*1024^2)

    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)
  }
}
