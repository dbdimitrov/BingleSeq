#' Bulk Functional Annotation Tab UI
#'
#' @export
#' @return None
sc_faUI <- function(id) {
  ns <- NS(id)
  
  tagList(# Sidebar panel for inputs ----
          sidebarPanel(tabsetPanel(
            id = ns("faSideTabSet"),
            
            tabPanel(
              value = ("faGetDorothea"),
              title = "TF activities",
              
              h4("Estimate TF Activities with DoRothEA"),
              
              tags$hr(),
              
              checkboxGroupInput(ns("faCheckBox"),
                                 "DoRothEA Confidence levels",
                                 c("A","B","C","D","E"),
                                 selected = c("A","B","C")),
              
              numericInput(
                ns("faTermRegulonValue"),
                label = "Minimum Regulon Size",
                value = 5,
                min = 3
              ),
              fluidRow(column(3, verbatimTextOutput(
                ns("faTermPvalue")
              ))),
              
              selectInput(
                ns("faOrganism"),
                label = "Select Organism",
                choices = list(
                  "Homo sapiens" = "hh",
                  "Mus musculus" = "mm"
                )
              ),
              
              numericInput(
                ns("faTFNumValue"),
                label = "Number of Top TFs to be plotted",
                value = 20,
                min = 5
              ),
              
              fluidRow(column(3, verbatimTextOutput(
                ns("faTFNumValue")))),
              
              actionButton(ns("faTFButton"),
                           label = "Get TF activities"),
              tags$hr()
            )
            ,
            
            
            tabPanel(
              value = ("faGetPROGENy"),
              title = "Pathway Activities",
              
              h4("Get Pathway Activities with PROGENy"),
              
              checkboxInput(ns("faProGeneType"),
                            label = ("Gene IDs are Symbols"),
                            value = TRUE),
              
              conditionalPanel(condition = "!input.faProGeneType",
                               ns = ns,
                               textInput(ns("faProGeneText"),
                                         "Enter Gene Type:",
                                         "ENSEMBL")
              ),
              
              tags$hr(),
              
              selectInput(
                ns("faProOrganism"),
                label = "Select Organism",
                choices = list(
                  "Homo sapiens" = "Human",
                  "Mus musculus" = "Mouse"
                )
              ),
              
              numericInput(
                ns("faProTopGeneNum"),
                label = "Number of footprint genes per pathway",
                value = 500,
                min = 10
              ),
              
              tags$hr(),
              
              actionButton(ns("faProgenyButton"),
                           label = "Get Pathway Activities per Sample"),
              
              
              tags$hr(),
              
            )
          )
          ),
          
          # Main panel for displaying outputs ----
          mainPanel(
            htmlOutput(ns("faInfo")),
            textOutput(ns("faTableText")),
            DT::dataTableOutput(ns("faTable")),
            
            conditionalPanel(condition = "input.faTFButton > 0",
                             ns = ns,
                             downloadButton(ns("faDownload"), "Download Table")
            ),
            tags$hr(),
            textOutput(ns("faPlotText")),
            plotOutput(ns("faPlot"), width = "800px", height = "700px")
          )
  )
}


#' Single cell FootPrint Tab Server
#'
#' @param finalData Seurat object with dim. red. and clustering results
#' @return None
sc_fa_Server <- function(input, output, session, finalData) {
  fa <- reactiveValues()
  
  output$faInfo <- renderUI({
    if(input$faTFButton == 0 && input$faProgenyButton == 0){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px'>
        This tab enables the use of footprint analysis tools
         DoRothEA and PROGENy which  infer the activity of TFs and pathways,
         respectively. Footprint-based strategies such as the aforementioned
         packages infer TF/pathway activityfrom the expression of molecules 
         considered to be donstream of a given pathway/TF.
         <br> <br>
         DoRothEA is a gene set resource containing signed TF-target 
         interactions that can be coupled with different statistical methods to
         estimate TF activity. In BingleSeq, DoRothEA is coupled to the
         statistical method VIPER.
         <br>
         'DoRothEA Confidence levels' are based on the supporting evidence for
         the TF-target interactions with A being the highest level.
         <br>
         'Minimum Regulon Size' refers to the minimum number of TF-target
         interactions for a given TF.
         <br> <br>
        PROGENy estimates the activity of 14 signalling pathways from gene
        expression using a linear model.
        PROGENy is based on downstream gene signatures observed to be 
        consistently deregulated in pertrubation experiments.
        
        'PROGENy footprint genes' refers to the number of most responsive genes
        per pathway. This number should be increased in cases with low gene
        coverage such as RNA-Seq (for more information refer to the papers cited
        in BingleSeq's manuscript).
         
         <br> <br>
         
         Also, please ensure that the approprate gene nomelcature is selected.

        </div>")
    }else {
      HTML("")
    }
  })
  
  observeEvent(input$faTFButton,{
    
    show_waiter(tagList(spin_folding_cube(), h2("Loading...")))
    
    fa$regulons <- fetch_regulons(input$faOrganism, input$faCheckBox)
    
    fa$dorothea_tfs <- sc_dorothea_viper(finalData$finalData,
                                         input$faTFNumValue,
                                         fa$regulons)
    
    waiter_hide()
    
    output$faTable <-
      DT::renderDataTable(fa$dorothea_tfs %>%
                            as.data.frame() %>%
                            t())
    
    output$faPlot <- renderPlot({
      sc_dorothea_plot(fa$dorothea_tfs)
    })
    
    output$faTableText <- renderText({
      "Normalized Enrichment Scores (NES) for each TF per cell cluster"
    })
    
    output$faPlotText <- renderText({
      "Heatmap with TF activity per cluster"
    })
  })
  
  # Progeny ----
  
  observeEvent(input$faProgenyButton,{
    
    show_waiter(tagList(spin_folding_cube(), h2("Loading...")))
    
    fa$pathways <- sc_progeny(finalData$finalData,
                              input$faProOrganism,
                              input$faProTopGeneNum,
                              1)
    
    waiter_hide()
    
    output$faTable <-
      DT::renderDataTable(fa$pathways %>%
                            as.data.frame() %>%
                            t())
    
    output$faPlot <- renderPlot({
      sc_progeny_plot(fa$pathways)
    })
    
    output$faTableText <- renderText({
      "Pathway activities represented by
      Normalised Enrichment Scores (NES) per Cell Cluster"
    })
    
    output$faPlotText <- renderText({
      "Normalised Enrichment Scores (NES) for each pathway per Cell Cluster"
    })
  })
}




#' Dorothea-Viper TF Activity estimation
#'
#' @param s_object Seurat object with clustering results
#' @param tf_num Minimum Regulon size
#' @param regulon dorothea regulon (human or mouse)
#' 
#' @keywords internal
#' @return Returns a table with TF activity scores ready for display
sc_dorothea_viper <- function(s_object, tf_num, regulon) {
  
  ## Compute Viper Scores
  s_object <- run_viper(s_object, regulon,
                        options = list(method = "scale", minsize = 5,
                                       eset.filter = FALSE, cores = 4,
                                       verbose = FALSE))
  
  DefaultAssay(object = s_object) <- "dorothea"
  s_object <- ScaleData(s_object)
  
  viper_scores_df <- 
    GetAssayData(s_object, slot = "scale.data", assay = "dorothea") %>%
    data.frame() %>%
    t()
  
  ## Create data frame with cells and clusters
  CellsClusters <- data.frame(cell = names(Idents(s_object)),
                              cell_type = as.character(Idents(s_object)),
                              stringsAsFactors = FALSE)  %>%
    mutate(cell = gsub('\\.', '-', cell))
  
  ## Create data frame with the Viper score per cell and its clusters
  viper_scores_clusters <- viper_scores_df  %>%
    data.frame() %>%
    rownames_to_column("cell") %>%
    mutate(cell = gsub('\\.', '-', cell))%>%
    gather(tf, activity, -cell) %>%
    inner_join(CellsClusters, by="cell")
  
  ## Summarize Viper scores by cell population
  summarized_viper_scores <- viper_scores_clusters %>%
    group_by(tf, cell_type) %>%
    summarise(avg = mean(activity),
              std = sd(activity))
  
  # Select the most variable TFs across clusters according to our scores.
  # (n Most variable TFs * number of clusters)
  tf_num.multiplied = tf_num * nlevels(s_object)
  
  highly_variable_tfs <- summarized_viper_scores %>%
    group_by(tf) %>%
    mutate(var = var(avg))  %>%
    ungroup() %>%
    top_n(tf_num.multiplied, var) %>%
    distinct(tf)
  
  # Table to be returned and plotted
  summarized_tfs <- summarized_viper_scores %>%
    semi_join(highly_variable_tfs, by = "tf") %>%
    dplyr::select(-std) %>%
    spread(tf, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  return(summarized_tfs)
}


#' Dorothea-Viper TF Activity Plot
#'
#' @param summarized_viper_scores_df a summarized df with TF activities
#' 
#' @keywords internal
#' @return Returns a heatmap with TF activities
sc_dorothea_plot <- function(summarized_viper_scores_df){
  palette_length = 100
  my_color = colorRampPalette(c("Darkblue", "white","red"))(palette_length)
  
  my_breaks <- c(seq(min(summarized_viper_scores_df), 0,
                     length.out=ceiling(palette_length/2) + 1),
                 seq(max(summarized_viper_scores_df)/palette_length,
                     max(summarized_viper_scores_df),
                     length.out=floor(palette_length/2)))
  
  dorothea_hmap <- pheatmap(t(summarized_viper_scores_df),fontsize=18,
                            fontsize_row = 14,
                            color=my_color, breaks = my_breaks,
                            main = "", angle_col = 45,
                            treeheight_col = 0,  border_color = NA)
  return(dorothea_hmap)
}



#' Progeny Pathway Activity estimation
#'
#' @param s_object Seurat object with clustering results
#' @param organism Human or Mouse Model Organism
#' @param top No. of most variable genes per pathway 
#' @param perm No. of permutations to be performed
#' 
#' @keywords internal
#' @return Returns a table with Pathway activity scores ready for display
sc_progeny <- function(s_object, organism, top, perm) {
  CellsClusters <- data.frame(Cell = names(Idents(s_object)), 
                              CellType = as.character(Idents(s_object)),
                              stringsAsFactors = FALSE)
  
  # Compute Progeny activity scores and assign to an assay 
  s_object <- progeny(s_object,
                      scale=FALSE,
                      organism=organism, top=top,
                      perm=perm, 
                      return_assay = TRUE)
  
  # Scale pathway activity scores. 
  s_object <- Seurat::ScaleData(s_object, assay = "progeny") 
  
  # Transform Progeny scores into a data frame
  progeny_scores_df <- 
    as.data.frame(t(GetAssayData(s_object, slot = "scale.data", 
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 
  
  # Match Progeny scores with cell clusters
  progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
  
  # Summarize Progeny scores by cell population
  summarized_progeny_scores <- progeny_scores_df %>% 
    group_by(Pathway, CellType) %>%
    summarise(avg = mean(Activity), std = sd(Activity))
  
  # Prepare the data for display
  progeny_scores_df <- summarized_progeny_scores %>%
    dplyr::select(-std) %>%   
    spread(Pathway, avg) %>%
    data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
  
  return(progeny_scores_df)
}


#' Progeny Pathway Activity Heatmap
#'
#' @param progeny_scores_df a summarized df with Pathway activities
#' 
#' @keywords internal
#' @return Returns a heatmap with TF activities
sc_progeny_plot <- function(progeny_scores_df){
  paletteLength = 100
  myColor = colorRampPalette(c("Darkblue", "white","red"))(paletteLength)
  
  progenyBreaks = c(seq(min(progeny_scores_df), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(progeny_scores_df)/paletteLength, 
                        max(progeny_scores_df), 
                        length.out=floor(paletteLength/2)))
  progeny_hmap = pheatmap(t(progeny_scores_df[,-1]), fontsize=18, 
                          fontsize_row = 14, 
                          color=myColor, breaks = progenyBreaks, 
                          main = "", angle_col = 45,
                          treeheight_col = 0,  border_color = NA)
  
  return(progeny_hmap)
}
