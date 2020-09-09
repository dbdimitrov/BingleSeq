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
                               h5("Available inputs"),
                               textOutput(ns("geneTypesText"))
                               
              ),
              
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
                  "Homo Sapiens" = "hh",
                  "Mus musculus" = "mm"
                )
              ),
              
              numericInput(
                ns("faTFNumValue"),
                label = "Number of Top TFs to be Plotted",
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
          )),
          
          # Main panel for displaying outputs ----
          mainPanel(
            tabsetPanel(
              tabPanel(
                value = ("faTable"),
                title = "Activities Table",
                
                htmlOutput(ns("faInfo")),
                
                DT::dataTableOutput(ns("faTable")),
                
                conditionalPanel(condition = "input.faTFButton > 0",
                                 ns = ns,
                                 downloadButton(ns("faDownload"), "Download Table")
                )
              ),
              
              tabPanel(
                value = ("faBarTab"),
                title = "Footprint Plot Tab",
                plotOutput(ns("faPlot"), width = "800px", height = "700px")
              )
            )
          ))
}


#' Bulk FootPrint Tab Server
#'
#' @param counts Unfiltered Count Table (Reactive Value)
#' @param de Differential Expression Results (Reactive Value)
#' @return None
bulk_faData <- function(input, output, session, counts, de) {
  fa <- reactiveValues()
  
  output$faInfo <- renderUI({
    if(input$faTFButton == 0){
      HTML("<div style='border:2px solid blue; font-size: 14px;
        padding-top: 8px; padding-bottom: 8px; border-radius: 10px'>
        info goes here
        </div>")
    }else {
        HTML("")
    }
  })
  
  output$geneTypesText <- renderPrint({
      print(AnnotationDbi::keytypes(org.Hs.eg.db))
  })
      
  observeEvent(input$faGeneType,{
      output$geneTypesText <- renderPrint({
        print(AnnotationDbi::keytypes(org.Hs.eg.db))
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
   
  })
  
  
  observeEvent(input$faSampleButton,{
    
    fa$tfs_sample <- fetch_tfs_per_sample(fa$data_de, input$faTermRegulonValue, fa$regulons, session)
    
    output$faTable <-
      DT::renderDataTable(as.data.frame(fa$tfs_sample), colnames = ("NES"))
    
    print(head(fa$tfs))
    print(head(fa$tfs_sample))
    
    output$faPlot <- renderPlot({
      plot_tfa_per_sample(fa$tfs, fa$tfs_sample, input$faTFNumValue)
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
    de_data$stat <- de_data[[3]]
    
    de_data_matrix <- de_data %>% 
      dplyr::select(gene, stat) %>% 
      dplyr::filter(!is.na(stat)) %>% 
      column_to_rownames(var = "gene") %>%
      as.matrix()
    
    out <- dorothea::run_viper(de_data_matrix, regulons,
                               options =  list(minsize = reg_size, eset.filter = FALSE, 
                                               cores = 1, verbose = FALSE, nes = TRUE))
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
           theme(axis.title = element_text(face = "bold", size = 12),
                 axis.text.x = 
                   element_text(angle = 45, hjust = 1, size =10, face= "bold"),
                 axis.text.y = element_text(size = 10, face= "bold"),
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
                          options =  list(minsize = reg_size, eset.filter = FALSE, 
                                          cores = 1, verbose = FALSE, method = c("scale")))
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
                            fontsize=14, fontsize_row = 12, fontsize_col = 12, 
                            color=myColor, breaks = dorotheaBreaks,
                            main = "", angle_col = 45,
                            treeheight_col = 0,  border_color = NA)
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
