#' Initialize Server
#'
#' @export
#' @return None
server <- function(input, output, session) {
  globalRV <- reactiveValues(app = NULL,
                             saveDir = NULL)


  observeEvent(input$saveFilesDirButton, {
    volumes <- getVolumes()

    shinyDirChoose(input,
                   'saveFilesDirButton',
                   roots = volumes,
                   session = session)

    savePath <- reactive({
      return(print(parseDirPath(volumes, input$saveFilesDirButton)))
    })


    # show meta data table
    if (!is.null(savePath)) {
      tryCatch({
        req(nchar(savePath()) > 0)

        globalRV$saveDir <- (data.dir = savePath())


      })

      # Create directories to save files
      dir.create(file.path(globalRV$saveDir, "figures"), showWarnings = FALSE)
      dir.create(file.path(globalRV$saveDir, "output"), showWarnings = FALSE)

    }
  })


  # Launch the app
  observeEvent(input$launch_app, {
    # if (!is.null(globalRV$saveDir)) {
      globalRV$app <- input$chooseApp

      removeTab("mainPage", "startApp", session = session)
    # } else{
    #
    #   sendSweetAlert(
    #     session = session,
    #     title = "No Output Directory",
    #     text = "Please select a directory using the 'Select Directory' button",
    #     type = "info"
    #   )
    # }
  })


  observe({


      if (!is.null(globalRV$app) && globalRV$app == 1) {
      ### Single-Cell RNA-Seq Pipeline App ----

        # Used to restrict tab creation
        tabs <- reactiveValues(
          counts = FALSE,
          filt = FALSE,
          norm = FALSE,
          clust = FALSE,
          markers = FALSE
        )


        # Append Single Cell Load Tab UI and swap ----
        appendTab(inputId = "mainPage",
                  tabPanel(
                    title = "Load Data",
                    sc_loadDataUI("loadTab"),
                    value = "loadTab"
                  ))

        updateTabsetPanel(session, 'mainPage', 'loadTab')


        # Load Data ----
        counts <- callModule(sc_loadData, "loadTab")

        observe({
          if (!is.null(counts$countTable) && !tabs$counts) {
            appendTab(
              inputId = "mainPage",
              tabPanel(
                id = "qcPage",
                value = "qcPage",
                title = "Quality Control",
                sc_qcUI("filterTab")
              )
            )
            tabs$counts = TRUE
          }
        })

        # Filter Data ----
        filteredData <- callModule(sc_qc, "filterTab", counts)

        observe({
          if (!is.null(filteredData$filteredData) && !tabs$filt) {
            appendTab(inputId = "mainPage",
                      tabPanel(title = "Normalization", sc_normUI("normTab")))

            tabs$filt = TRUE
          }
        })


        # Normalize Data ----
        normalizedData <-
          callModule(sc_norm, "normTab", filteredData)

        observe({
          if (!is.null(normalizedData$normalizedData) && !tabs$norm) {
            appendTab(inputId = "mainPage",
                      tabPanel(title = "Clustering", sc_clustUI("clustTab")))

            tabs$norm = TRUE
          }
        })


        # Cluster Data ----
        finalData <-
          callModule(sc_clust, "clustTab", normalizedData)

        observe({
          if (!is.null(finalData$finalData) && !tabs$clust) {
            appendTab(inputId = "mainPage",
                      tabPanel(title = "Differential Expression", sc_deUI("deTab")))
            appendTab(inputId = "mainPage",
                      tabPanel(title = "Compare DE Methods", sc_compUI("compTab")))

            tabs$clust = TRUE
          }
        })

        observe({
          if (!is.null(finalData$finalData)) {

            updateNumericInput(
              session = session,
              inputId = "deTab-dgeClustInput",
              label = "Select Cluster of Interest",
              value = min(as.numeric(
                levels(finalData$finalData@active.ident)
              )),
              min = min(as.numeric(
                levels(finalData$finalData@active.ident)
              )),
              max = max(as.numeric(
                levels(finalData$finalData@active.ident)
              ))
            )
          }
        })



        # Differential Expression Markers ----
        markers <- callModule(sc_de, "deTab", finalData)

        observe({
          if (!is.null(markers$markers) && !tabs$markers) {
            appendTab(inputId = "mainPage",
                      tabPanel(title = "Functional Annotation", sc_goUI("goTab")))

            tabs$markers = TRUE
          }

          if (!is.null(finalData$finalData)) {
            # Provide min and max clust No for plots
            output$goClustNoInputUI <- renderUI({
              numericInput(("goTab-goClustNoInput"),
                           label = "Select Cluster of Interest",
                           value = min(as.numeric(
                             levels(finalData$finalData@active.ident)
                           )),
                           min = min(as.numeric(
                             levels(finalData$finalData@active.ident)
                           )),
                           max = max(as.numeric(
                             levels(finalData$finalData@active.ident)
                           ))
              )

            })
          }
        })


        # Functional Annotation ----
        callModule(sc_go, "goTab", markers, counts)



        # Compare DE Methods ----
        callModule(sc_comp, "compTab", finalData)


      } else if (!is.null(globalRV$app) && globalRV$app == 2) {
      ### Bulk RNA-Seq Pipeline App ----


        tabs <- reactiveValues(# Used to restrict tab creation
          counts = FALSE,
          filt = FALSE,
          de = FALSE)


        plotChoices = list(
          "PCA plot" = "pca",
          "Scree plot" = "scree",
          "Barchart" = "bar",
          "Volcano Plot" = "volcano",
          "MA Plot" = "MA",
          "Heatmap" = "heat"
        )

        vennChoice = list("Venn Diagram" = "venn")



        # Append Bulk Load Tab and Swap ----
        appendTab(inputId = "mainPage",
                  tabPanel(
                    title = "Load Data",
                    bulk_loadDataUI("loadTab"),
                    value = "loadTab"
                  ))

        updateTabsetPanel(session, 'mainPage', 'loadTab')



        # Load Data ----
        counts <- callModule(bulk_loadData, "loadTab")

        observe({
          if (!is.null(counts$countTable) && !tabs$counts) {
            appendTab(
              inputId = "mainPage",
              tabPanel(
                id = "filterTab",
                value = "filterTab",
                title = "Quality Control",
                bulk_filterDataUI("filterTab")
              )
            )

            tabs$counts = TRUE
          }

          if (!is.null(counts$countTable)) {
            # Filter Tab Table prefilter (to be changed) ----
            output$"filterTab-prefilterTable" <-
              DT::renderDataTable(DT::datatable(
                generateSummary(counts$countTable, session),
                options = list(paging = FALSE, searching = FALSE),
                rownames = FALSE
              ))
          }
        })


        # Filter Data ----
        filt <- callModule(bulk_filterData, "filterTab", counts)

        observe({
          if (!is.null(filt$filteredCounts) && !tabs$filt && filt$correctFormat) {
            appendTab(
              inputId = "mainPage",
              tabPanel(
                id = "deTab",
                value = "deTab",
                title = "Differential Expression",
                bulk_deDataUI("deTab")
              )
            )

            tabs$filt = TRUE
          }
        })


        # Differential Expression Data -----
        de <- callModule(bulk_deData, "deTab", filt)

        observe({
          if (!is.null(de$deTable) && !tabs$de) {
            appendTab(
              inputId = "mainPage",
              tabPanel(
                id = "plotTab",
                value = "plotTab",
                title = "Visualize data",
                bulk_plotDataUI("plotTab")
              )
            )

            appendTab(
              inputId = "mainPage",
              tabPanel(
                id = "goTab",
                value = "goTab",
                title = "Functional Annotation",
                bulk_goDataUI("goTab")
              )
            )

            appendTab(
              inputId = "mainPage",
              tabPanel(
                id = "comTab",
                value = "compTab",
                title = "DE Package Comparison",
                bulk_compDataUI("compTab")
              )
            )

            tabs$de = TRUE
          }
        })

        # Dynamically show/hide condition options
        observeEvent(req(input$mainPage) == "plotTab", {
          if (!is.null(de$deTable) && de$conditionNo > 2) {
            shinyjs::show("plotTab-selectConditionMA")
            shinyjs::show("plotTab-selectConditionTHM")
            shinyjs::show("plotTab-goGetCondition")

            updateSelectInput(session,
                              "plotTab-selectPlotCombo",
                              choices = append(plotChoices, vennChoice))

          } else if (!is.null(de$deTable) && de$conditionNo == 2) {
            updateSelectInput(session, "plotTab-selectPlotCombo", choices = plotChoices)

            print("resets")

            updateSelectInput(session, "plotTab-selectConditionMA", selected = 1)
            updateSelectInput(session, "plotTab-selectConditionTHM", selected = 1)
            updateSelectInput(session, "goTab-goGetCondition", selected = 1)

            shinyjs::hide("plotTab-selectConditionMA")
            shinyjs::hide("plotTab-selectConditionTHM")
            shinyjs::hide("goTab-goGetCondition")

          }
        })


        # Plot Data ------
        callModule(bulk_plotData, "plotTab", de)

        # Functional Annotation -----
        callModule(bulk_goData, "goTab", counts , de)

        # Compare Data ------
        callModule(bulk_compData, "compTab", filt, de)
      }
  })
}
