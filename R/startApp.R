#' Launch the application
#'
#' #' Launches the app
#'
#' @import shiny
#' @import shinyjs
#' @import shinyFiles
#' @import shinyWidgets
#' @import V8
#' @import DT
#' @import grDevices
#' @import ggplot2
#' @import VennDiagram
#' @import factoextra
#' @import fastcluster
#' @import reshape2
#' @import gridExtra
#' @import ggrepel
#' @import DESeq2
#' @import edgeR
#' @import limma
#' @import ALDEx2
#' @import MAST
#' @import monocle
#' @import Seurat
#' @import scran
#' @import SC3
#' @import goseq
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import GO.db
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

    library(V8) # enables js lines to be read intext
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

    library(DESeq2)
    library(edgeR)
    library(limma)
    library(ALDEx2)
    library("MAST")

    # Load required SC packages
    library(monocle)
    library(Seurat)
    library('scran')
    library('SC3')

    library(goseq)
    library(org.Mm.eg.db)
    library(org.Hs.eg.db)
    library(GO.db)

    app <- shinyApp(ui = shinyUI(ui),
                    server = shinyServer(server))
    runApp(app)
  }
}
