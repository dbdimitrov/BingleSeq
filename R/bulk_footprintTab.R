#' Bulk Functional Annotation Tab UI
#'
#' @export
#' @return None
bulk_faDataUI <- function(id) {
  ns <- NS(id)
  
  tagList(# Sidebar panel for inputs ----
          sidebarPanel(tabsetPanel(
            id = ns("faSideTabSet"),
            
            tabPanel(
              value = ("faGetDorothea"),
              title = "TF activities",
              
              h4("Get Transcription Factor Activities with DoRothEA"),
              
              checkboxInput(ns("faGeneType"), label = ("Gene IDs are Symbols"),
                            value = TRUE),
              
              conditionalPanel(condition = "!input.faGeneType",
                               ns = ns,
                               textInput(ns("faTypeInput"),
                                         "Enter Gene Type:",
                                         "ENSEMBL"),
                               tags$b("Available inputs:"),
                               textOutput(ns("geneTypesText"))
                               
              ),
              
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
                value = 25,
                min = 5
              ),
              fluidRow(column(3, verbatimTextOutput(
                ns("faTFNumValue")))),
              
              actionButton(ns("faTFButton"),
                           label = "Get TF activities"),
              conditionalPanel(condition = "input.faTFButton>0",
                               ns = ns,
                               actionButton(ns("faSampleButton"),
                                            label = "Get TF activities per Sample")
              ),
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
                value = 100,
                min = 10
              ),
              
              tags$hr(),
              
              actionButton(ns("faProSampleButton"),
                           label = "Get Pathway Activities per Sample"),
              
              tags$hr(),
              
              actionButton(ns("faProConditionButton"),
                           label = "Get Pathway Activities"),
              
              
              tags$hr(),
              
              conditionalPanel(condition = "input.faProConditionButton > 0",
                               ns = ns,

                               tags$b("Available pathways:"),
                               tags$i("Androgen, EGFR, Estrogen, Hypoxia,
                                         JAK-STAT, MAPK, NFkB, p53, PI3K, TGFb,
                                         TNFa, Trail, VEGF, WNT"),
                               textInput(ns("faProPathwayText"),
                                         "Enter Pathway:",
                                         "JAK-STAT"),
                               actionButton(ns("faProPathwayButton"),
                                            label = "Get Top Genes per Pathway")
                               
                               
              ),
              
            )
          )
        ),
          
          # Main panel for displaying outputs ----
          mainPanel(
              htmlOutput(ns("faInfo")),
              textOutput(ns("faPlotText")),
              DT::dataTableOutput(ns("faTable")),
              
              conditionalPanel(condition = "input.faTFButton > 0",
                               ns = ns,
                               downloadButton(ns("faDownload"), "Download Table")
              ),
              tags$hr(),
              textOutput(ns("faTableText")),
              plotOutput(ns("faPlot"), width = "800px", height = "700px")
          )
        )
}


#' Bulk FootPrint Tab Server
#'
#' @param counts Unfiltered Count Table (Reactive Value)
#' @param de Differential Expression Results (Reactive Value)
#' @return None
bulk_faData <- function(input, output, session, counts, de) {
  fa <- reactiveValues()
  
  output$faInfo <- renderUI({
    if(input$faTFButton == 0 && input$faProSampleButton == 0
       && input$faProConditionButton == 0){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px'>
        This tab enables the use of footprint analysis tools
        DoRothEA and PROGENy which are tools used infer the activity of TFs and
        pathways, respectively. Footprint-based strategies such as the
        aforementioned packages infer TF/pathway activity from the expression
        of molecules considered to be downstream of a given pathway/TF.
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
  
  
  output$geneTypesText <- renderPrint({
    "ALIAS, ENSEMBL, ENTREZID, ENZYME, SYMBOL, UNIGENE"
  })
      
  observeEvent(input$faGeneType,{
      output$geneTypesText <- renderPrint({
        "ALIAS, ENSEMBL, ENTREZID, ENZYME, SYMBOL, UNIGENE"
      })
  })
    
  observeEvent(input$faTFButton,{
    
   fa$data_de <- tibble::rownames_to_column(de$merged, "X")
   if(!input$faGeneType){
     fa$data_de <- convert_gene_type(fa$data_de, input$faTypeInput, "SYMBOL")
   }
   fa$regulons <- fetch_regulons(input$faOrganism, input$faCheckBox)
   fa$tfs <- fetch_tf_activities(fa$data_de, input$faTermRegulonValue, fa$regulons, session)
   
   output$faTable <-
     DT::renderDataTable(as.data.frame(fa$tfs), colnames = ("NES"))
   
   output$faPlot <- renderPlot({
     plot_top_tfs(fa$tfs, input$faTFNumValue)
   })
   
   output$faTableText <- renderText({
     "Table with TFs and corresponding NES activity estimate"
   })
   
   output$faPlotText <- renderText({
     "Top TFs according to normalized enrichment scores (NES)"
   })
   
   
  })
  
  
  observeEvent(input$faSampleButton,{
    
    fa$tfs_sample <- fetch_tfs_per_sample(fa$data_de,
                                          input$faTermRegulonValue,
                                          fa$regulons,
                                          session)
    
    output$faTable <-
      DT::renderDataTable(as.data.frame(fa$tfs_sample))
    
    
    output$faPlot <- renderPlot({
      plot_tfa_per_sample(fa$tfs, fa$tfs_sample, input$faTFNumValue)
    })
    
    output$faTableText <- renderText({
      "Table with TFs and corresponding NES per Sample"
    })
    
    output$faPlotText <- renderText({
      "Top TFs with corresponding NES per Sample"
    })
  })
  
  
  # Progeny ----
  
  observeEvent(input$faProSampleButton,{
    
    fa$data_de <- tibble::rownames_to_column(de$merged, "X")
    if(!input$faProGeneType){
      fa$data_de <- convert_gene_type(fa$data_de, input$faTypeInput, "SYMBOL")
    }
    fa$progSample<- get_pathway_activity_per_sample(fa$data_de,
                                                    as.character(input$faProOrganism),
                                                    as.numeric(input$faProTopGeneNum),
                                                    session)

    output$faTable <-
      DT::renderDataTable(as.data.frame(fa$progSample[[1]]))

    output$faPlot <- renderPlot({
      grid.draw(fa$progSample[[2]])
    })
    
    output$faTableText <- renderText({
      "Table with Pathway activities represented by NES per Sample"
    })
    
    output$faPlotText <- renderText({
      "Pathway activities per Sample Heatmap"
    })
  })
  
  
  observeEvent(input$faProConditionButton,{
    
    fa$data_de <- tibble::rownames_to_column(de$merged, "X")
    if(!input$faProGeneType){
      fa$data_de <- convert_gene_type(fa$data_de, input$faTypeInput, "SYMBOL")
    }
    fa$proCondition <- get_pathway_activity(fa$data_de,
                                             as.character(input$faProOrganism),
                                             as.numeric(input$faProTopGeneNum),
                                             session)
    
    output$faTable <-
      DT::renderDataTable(as.data.frame(fa$proCondition[[1]]))
    
    output$faPlot <- renderPlot({
      grid.draw(fa$proCondition[[2]])
    })
    
    output$faTableText <- renderText({
      "Table with Pathway activities represented by NES per Sample"
    })
    
    output$faPlotText <- renderText({
      "Normalised Enrichment Scores (NES) for each pathway"
    })
  })
  
  
  observeEvent(input$faProPathwayButton,{
    
    if(req(fa$proCondition[[3]])){
      fa$scatter <- get_progeny_scatter(fa$proCondition[[3]],
                                        as.character(input$faProPathwayText),
                                        as.character(input$faProOrganism),
                                        as.numeric(input$faProTopGeneNum),
                                        fa$de_data
      )
    } 
    plot(fa$scatter)
    output$faPlot <- renderPlot({
      fa$scatter
    })
    
    output$faPlotText <- renderText({
      "Scatter plot with n most responsive genes for the provided pathway"
    })
  })
}



#' Transcription Factor activity with DoRothEA
#'
#' @param de_data Table with de results and symbols for gene type
#' @param reg_size minimun size of regulons
#' @param regulons Dorothea regulons
#' @param session currently running r session
#' @export
#' @return Returns a data frame with TF activities estimated from TF footprints
fetch_tf_activities <- function(de_data, reg_size, regulons, session) {
  out <- tryCatch({
    de_data$gene <- de_data$X
    gene_counts <- de_data[,6:ncol(de_data)]
    de_data$stat <- de_data[[3]] # reassign package-specific stat 
    
    de_data_matrix <- de_data %>% 
      dplyr::select(gene, stat) %>% 
      dplyr::filter(!is.na(stat)) %>% 
      column_to_rownames(var = "gene") %>%
      as.matrix()
    
    out <- dorothea::run_viper(de_data_matrix, regulons,
                               options =  list(minsize = reg_size,
                                               eset.filter = FALSE, 
                                               cores = 1, verbose = FALSE,
                                               nes = TRUE))
  },
  error = function(cond) {
    sendSweetAlert(
      session = session,
      title = "Data Mismatch",
      text = "Please ensure that Gene IDs are of type Symbol, if select appropriate type to convert",
      type = "error"
    )
    return()
  })
  return(out)
}

#' Plot the top n NES for the TFs
#'
#' @param tf_activities a df with TF activities estimated from TF footprints
#' @param tf_num number of TFs with top absolute NES to be displayed
#' @export
#' @return 
plot_top_tfs <- function(tf_activities, tf_num) {
  
  tf_activities_top <- tf_activities %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "stat") %>%
    dplyr::top_n(tf_num, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))
  
  
  return(ggplot(tf_activities_top,aes(x = reorder(GeneID, NES), y = NES)) + 
           geom_bar(aes(fill = NES), stat = "identity") +
           scale_fill_gradient2(low = "darkblue", high = "indianred", 
                                mid = "whitesmoke", midpoint = 0) + 
           theme_minimal() +
           theme(axis.title = element_text(face = "bold", size = 18),
                 axis.text.x = 
                   element_text(angle = 45, hjust = 1, size =16, face= "bold"),
                 axis.text.y = element_text(size = 16, face= "bold"),
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) +
           xlab("Transcription Factors"))
}




#' TF activity per sample Sample DoRothEA
#'
#' @param de_data Table with de results and symbols for gene type
#' @param reg_size minimum size of regulons
#' @param regulons Regulons obtained from dorothea
#' @export
#' @return Returns a data frame with TF activities per Sample
fetch_tfs_per_sample <- function(de_data, reg_size, regulons, session) {
  out <- tryCatch({
    normalised_counts <- de_data[,6:ncol(de_data)]
    normalised_counts$gene <- de_data$X
    
    normalised_counts_matrix <- normalised_counts %>% 
      dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
      tibble::column_to_rownames(var = "gene") %>% 
      as.matrix()
    
    out <- 
      dorothea::run_viper(normalised_counts_matrix, regulons,
                          options =  list(minsize = reg_size,
                                          eset.filter = FALSE, 
                                          cores = 1, verbose = FALSE,
                                          method = c("scale")))
  },
  error = function(cond) {
    sendSweetAlert(
      session = session,
      title = "TF activities not found",
      text = "Please ensure that overall TF activity is estimated",
      type = "error"
    )
    return()
  })
  return(out)
}

#' Plot the top n TF activities per Sample in a heatmap
#'
#' @param tf_activities a df with TF activities estimated from TF footprints
#' @param tf_activities_counts a df with TF per Sample
#' @param tf_num number of TFs with top absolute NES to be displayed
#' @export
#' @return 
plot_tfa_per_sample <- function(tf_activities, tf_activities_counts, tf_num) {
  tf_activities_top <- tf_activities %>%
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::rename(NES = "stat") %>%
    dplyr::top_n(tf_num, wt = abs(NES)) %>%
    dplyr::arrange(NES) %>% 
    dplyr::mutate(GeneID = factor(GeneID))
  
  tf_activities_counts_filter <- tf_activities_counts %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "GeneID") %>%
    dplyr::filter(GeneID %in% tf_activities_top$GeneID) %>%
    column_to_rownames(var = "GeneID") %>%
    as.matrix()
  
  tf_activities_vector <- as.vector(tf_activities_counts_filter)
  
  paletteLength <- 100
  myColor <- 
    colorRampPalette(c("darkblue", "white","red"))(paletteLength)
  
  dorotheaBreaks <- c(seq(min(tf_activities_vector), 0, 
                          length.out=ceiling(paletteLength/2) + 1),
                      seq(max(tf_activities_vector)/paletteLength, 
                          max(tf_activities_vector), 
                          length.out=floor(paletteLength/2)))
  dorothea_hmap <- pheatmap(tf_activities_counts_filter,
                            fontsize=18, fontsize_row = 16, fontsize_col = 16, 
                            color=myColor, breaks = dorotheaBreaks,
                            main = "", angle_col = 45,
                            treeheight_col = 0,  border_color = NA, silent = T)
  
  return(dorothea_hmap)
}



#' Gene type converter function
#'
#' @param de_data Table with de results
#' @param input_type input gene type e.g. ENSEMBL
#' @param output_type output gene type e.g. Symbol
#' @export
#' @return Returns a df with gene names convert to symbols
convert_gene_type <- function(de_data, input_type, output_type) {
  require(AnnotationDbi)
  geneIDs1 <- AnnotationDbi::select(org.Hs.eg.db,
                                    keys=de_data$X,
                                    keytype = input_type, 
                                    columns = c(output_type, input_type))
  
  geneIDs1 <- subset(geneIDs1, (!duplicated(geneIDs1[[output_type]])))
  geneIDs1 <- subset(geneIDs1, !is.na(geneIDs1[[output_type]]))
  data_me <- merge(de_data, geneIDs1, by.x = "X", by.y = input_type)
  data_me$X <- data_me[[output_type]]
  data_me <- within(data_me, rm(list=sub("[.]test","",output_type)))
  return(data_me)
}

#' Fetch Dorothea Regulon
#'
#' @param org Human or Mouse
#' @param conf_list list_with_conf
#' @export
#' @return Returns Dorothea regulons
fetch_regulons <- function(org, conf_list) {
  if(org=="hh"){
    data(dorothea_hs, package = "dorothea")
    regulons = dorothea_hs %>%
      dplyr::filter(confidence %in% conf_list)
  }
  else if(org=="mm"){
    data(dorothea_mm, package = "dorothea")
    regulons = dorothea_mm %>%
      dplyr::filter(confidence %in% conf_list)
  }
  
  return(regulons)
}




#' Fetch Pathway Activity per Sample
#' 
#' @param de_data de results dataframe
#' @param organism Human or Mouse
#' @param top # the top n genes for generating the model matrix sorted by p
#' @export
#' @return a list with pathway activity counts and a heatmap 
get_pathway_activity_per_sample <- function(de_data, organism, top, session) {
  out <- tryCatch({
    # get and format normalized counts
    normalised_counts <- de_data[,6:ncol(de_data)]
    normalised_counts$gene <- de_data$X
    
    normalised_counts_matrix <- normalised_counts %>% 
      dplyr::mutate_if(~ any(is.na(.x)),~ if_else(is.na(.x),0,.x)) %>% 
      tibble::column_to_rownames(var = "gene") %>% 
      as.matrix()
    
    
    # Pathway activity with Progeny
    pathway_activity_counts <- progeny(normalised_counts_matrix,
                                       scale=TRUE, 
                                       organism=as.character(organism),
                                       top = top)
    activity_counts <- as.vector(pathway_activity_counts)
    
    # heatmap
    paletteLength <- 100
    myColor <- 
      colorRampPalette(c("darkblue", "whitesmoke","indianred"))(paletteLength)
    
    progenyBreaks <- c(seq(min(activity_counts), 0, 
                           length.out=ceiling(paletteLength/2) + 1),
                       seq(max(activity_counts)/paletteLength, 
                           max(activity_counts), 
                           length.out=floor(paletteLength/2)))
    
    progeny_hmap <- pheatmap(t(pathway_activity_counts),fontsize=18, 
                             fontsize_row = 16, fontsize_col = 16, 
                             color=myColor, breaks = progenyBreaks, 
                             main = "", angle_col = 45,
                             treeheight_col = 0,  border_color = NA)
    
    out <- list(pathway_activity_counts, progeny_hmap)
  },
  error = function(cond) {
    sendSweetAlert(
      session = session,
      title = "Data Mismatch",
      text = "Please ensure that Gene IDs are of type Symbol, if select appropriate type to convert",
      type = "error"
    )
    return()
  })
  return(out)
}


#' Fetch Pathway Activity
#' 
#' @param de_data de results dataframe
#' @param organism Human or Mouse
#' @param top num of top n genes for generating the model matrix sorted by p
#' @param session current R session
#' @export
#' @return a list with pathway activities, barplot, de_results matrix
get_pathway_activity <- function(de_data, organism, top, session) {
  out <- tryCatch({
    de_results <- de_data[,0:5]
    de_results$stat <- de_data[[3]] # reasign package-specific stat 
    
    de_results_matrix <- de_results %>% 
      dplyr::select(X, stat) %>% 
      dplyr::filter(!is.na(stat)) %>% 
      column_to_rownames(var = "X") %>%
      as.matrix()
    
    pathway_activity_zscore <- progeny(de_results_matrix, 
                                       scale=TRUE,
                                       organism=as.character(organism),
                                       top = top, perm = 10000,
                                       z_scores = TRUE) %>%
      t()
    colnames(pathway_activity_zscore) <- "NES"
    
    
    pathway_activity_zscore_df <- as.data.frame(pathway_activity_zscore) %>% 
      rownames_to_column(var = "Pathway") %>%
      dplyr::arrange(NES) %>%
      dplyr::mutate(Pathway = factor(Pathway))
    
    
    pathway_plot <- 
      ggplot(pathway_activity_zscore_df,aes(x = reorder(Pathway, NES), y = NES)) + 
      geom_bar(aes(fill = NES), stat = "identity") +
      scale_fill_gradient2(low = "darkblue", high = "indianred", 
                           mid = "whitesmoke", midpoint = 0) + 
      theme_minimal() +
      theme(axis.title = element_text(face = "bold", size = 18),
            axis.text.x = 
              element_text(angle = 45, hjust = 1, size = 16, face= "bold"),
            axis.text.y = element_text(size = 16, face= "bold"),
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank()) +
      xlab("Pathways")
    
    out <- list(pathway_activity_zscore_df, pathway_plot, de_results_matrix)
  },
  error = function(cond) {
    sendSweetAlert(
      session = session,
      title = "Data Mismatch",
      text = "Please ensure that Gene IDs are of type Symbol, if select appropriate type to convert",
      type = "error"
    )
    return()
  })
  return(out)
}




#' Fetch Progeny Scatter with top Genes
#' 
#' @param de_data_results de results matrix as returned by get_pathway_activity
#' @param x_pathway pathway to be plotted
#' @param organism Human or Mouse
#' @param top # of significant genes for each pathway
#' @param session current R session
#' 
#' @export
#' @return a scatter plot with the top most signif signature genes
get_progeny_scatter <- function(de_results_matrix, x_pathway,
                                organism, top, de_data) {
  
  prog_matrix <- getModel(as.character(organism), top=top) %>% 
    as.data.frame()  %>%
    tibble::rownames_to_column("GeneID")
  
  ttop_results <- de_results_matrix %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column("GeneID")
  
  # statistic name
  stat_x = as.character(colnames(de_data)[3])
  
  scat_plots <- progeny::progenyScatter(df = ttop_results, 
                                        weight_matrix = prog_matrix, 
                                        statName = stat_x, verbose = TRUE)
  
  return(scat_plots[[1]][[as.character(x_pathway)]])
}
