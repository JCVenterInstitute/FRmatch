library(shiny)
library(plotly)
library(FRmatch)
library(dplyr)
library(tibble)
library(SingleCellExperiment)
# library(shinycssloaders)

data("sce-layer1-15clusters.rda")
data("sce-layer1-topNodes.rda")

# dev.off() #FIGURE OUT WHY!!!???

#######################################################################################################
## some useful function

myfun.datasplit <- function(sce.query, sce.ref, seed=999, frac.ref=.5){
  set.seed(seed)
  ## subsampling
  all <- colData(sce.ref) %>% as.data.frame() %>% rownames_to_column() %>% as_tibble()
  sam1 <- all %>% group_by(cluster_membership) %>% sample_frac(frac.ref)
  sam2 <- dplyr::setdiff(all, sam1)

  newsce.ref <- sce.ref[,sam1$rowname] #reference
  newsce.query <- sce.query[,sam2$rowname] #query
  return(list("newsce.query"=newsce.query, "newsce.ref"=newsce.ref))
}

#######################################################################################################

# Define UI for dataset viewer app ----
ui <- fluidPage(

  # App title ----
  titlePanel("FRmatch Demo"),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      ## Input: Select datasets ----
      selectInput("querydata", "Choose a query dataset:",
                  choices = c("Layer1", "Layer1-topNodes")),
      selectInput("refdata", "Choose a reference dataset:",
                  choices = c("Layer1", "Layer1-topNodes")),
      # Input: Checkbox if to impute reference ----
      checkboxInput("imputation", "Imputation", FALSE),
      ## Include clarifying text ----
      helpText(p("Note: by checking the Imputation box,",
               "FRmatch will impute the dropout values for each marker gene",
               em("only"), "in the cluster that it marks and",
               em("only"), "in the reference dataset.")),

      # Horizontal line ----
      tags$hr(),

      ## Input: Specify random seed ----
      numericInput("seed", "Please set a random seed:", 100),

      ## Input: Specify spliting fraction ----
      sliderInput("splitFrac", "Data subsampling fraction:",
                  0.2, 0.8, 0.5, step=0.1),
      ## Include clarifying text ----
      helpText("Note: the above fraction of cells will be selected",
               "in proportion to the cluster sizes from the reference dataset,",
               "and the rest of the cells will be selected from the query dataset."),

      # Horizontal line ----
      tags$hr(),

      ## Input: actionButton() to defer the rendering of output ----
      actionButton("updateButton", "Run FRmatch",
                   class = "btn-primary")



    ), #close sidebarPanel


    # Main panel for displaying outputs ----
    mainPanel(

      tabsetPanel(
        tabPanel("Data Subsampling",
                 h4("Query Data"),
                 verbatimTextOutput("viewQuery"),
                 tableOutput("tableQuery"),
                 h4("Reference Data"),
                 verbatimTextOutput("viewRef"),
                 tableOutput("tableRef")),
        tabPanel("Dropouts",
                 h4("Check dropouts in the reference data"),
                 checkboxInput("afterImputation", "Show after imputation", FALSE),
                 plotOutput("dropouts", height="auto")),
        tabPanel("Results",
                 h4("Recommended matches"),
                 plotOutput("matches"),
                 sliderInput("sigLvl", "Significance level:",
                             0.01, 0.2, 0.05, step=0.01, width="400px"),
                 helpText("Note: by setting smaller significance level,",
                          "there will be more matches found."),
                 hr(),
                 h4("P-values"),
                 plotOutput("pvals"),
                 h4("FR statistics"),
                 plotOutput("FRstats")),
        tabPanel("MST",
                 h4("Minimum spanning tree"),
                 uiOutput("querySelection"),
                 uiOutput("refSelection"),
                 plotOutput("temp")
                 # uiOutput("plotMST")
                 )

      ) #close tabsetPanel

    ) #closde mainPanel

  ) #close sidebarLayout

) #close fliudPage

#######################################################################################################

# Define server logic to summarize and view selected dataset ----
server <- function(input, output, session) {

  # Return the requested dataset ----
  # Note that we use eventReactive() here, which depends on
  # input$update (the action button), so that the output is only
  # updated when the user clicks the button
  queryInput <- eventReactive(input$updateButton, {
    switch(input$querydata,
           "Layer1" = sce.layer1.15clusters,
           "Layer1-topNodes" = sce.layer1.topNodes)
  }, ignoreNULL = FALSE)
  refInput <- eventReactive(input$updateButton, {
    switch(input$refdata,
           "Layer1" = sce.layer1.15clusters,
           "Layer1-topNodes" = sce.layer1.topNodes)
  }, ignoreNULL = FALSE)

  ## subsampling data
  newData <- eventReactive(input$updateButton, {
    myfun.datasplit(queryInput(), refInput(), seed=input$seed, frac.ref=input$splitFrac)
  }, ignoreNULL = FALSE)
  ## subsampled reference data
  newDataRef <- eventReactive(input$updateButton, {
    newData()$newsce.ref
  }, ignoreNULL = FALSE)
  ## subsampled query data
  newDataQuery <- eventReactive(input$updateButton, {
    newData()$newsce.query
  }, ignoreNULL = FALSE)

  ## run FRmatch
  results <- eventReactive(input$updateButton, {
    FRmatch(newDataQuery(), newDataRef(), imputation=input$imputation)
  }, ignoreNULL = FALSE)

  ##------ tab: Data Subsampling ------##

  ## look at query data
  output$viewQuery <- renderPrint({
    newDataQuery()
  })
  output$tableQuery <- renderTable({
    tab.query <- table(colData(queryInput())$cluster_membership)
    tab.newquery <- table(colData(newDataQuery())$cluster_membership)
    table.query <- cbind(tab.query, tab.newquery) %>% data.frame() %>% rownames_to_column()
    colnames(table.query) <- c("Cluster", "Size", "Selected")
    table.query
  }, rownames=TRUE)

  ## look at ref data
  output$viewRef <- renderPrint({
    newDataRef()
  })
  output$tableRef <- renderTable({
    tab.ref <- table(colData(refInput())$cluster_membership)
    tab.newref <- table(colData(newDataRef())$cluster_membership)
    table.ref <- cbind(tab.ref, tab.newref) %>% data.frame() %>% rownames_to_column()
    colnames(table.ref) <- c("Cluster", "Size", "Selected")
    table.ref
  }, rownames=TRUE)

  ##------ tab: Dropouts ------##

  ## dropout plot
  output$dropouts <- renderPlot({
    if(input$afterImputation){
      newDataRefImputation <- impute_dropout(newDataRef())
      check_dropout(newDataRefImputation, return.value=FALSE, plot.dropout=TRUE)
    }
    else check_dropout(newDataRef(), return.value=FALSE, plot.dropout=TRUE)
  }, height = function() {
    session$clientData$output_dropouts_width
  })

  ##------ tab: Results ------##

  ## matches plot
  output$matches <- renderPlot({
    plot_FRmatch(results(), sig.level=input$sigLvl)
  })

  ## p-values plot
  output$pvals <- renderPlot({
    plot_FRmatch(results(), type="pmat")
  })

  ## FR stats plot
  output$FRstats <- renderPlot({
    plot_FRmatch(results(), type="statmat")
  })

  ##------ tab: MST ------##
  output$querySelection <- renderUI({
    selectInput("queryCluster", "Select query cluster:", choices = newDataQuery()@metadata$cluster_order)
  })
  output$refSelection <- renderUI({
    selectInput("refCluster", "Select reference cluster:", choices = newDataRef()@metadata$cluster_order,
                multiple = TRUE)
  })
  ## This is the function to break the whole data into different blocks for each plot
  plotInput <- reactive({
    markergenes <- newDataRef()@metadata$cluster_marker_info$markerGene
    query.cluster <- input$queryCluster
    ind.query <- which(colData(newDataQuery())$cluster_membership==query.cluster)
    samp.query <- unname(logcounts(newDataQuery())[markergenes, ind.query])
    return(list("markergenes"=markergenes, "samp.query"=samp.query))
  })
  ## plots
  # output$plotMST <- renderUI({
  #   plot_output_list <- lapply(1:length(input$refCluster), function(i) {
  #     plotname <- paste("plot", i, sep="")
  #     plotOutput(plotname, width="400px", height="400px")
  #   })
  #   do.call(tagList, plot_output_list)
  # })
  # observe({
  #   lapply(1:length(input$refCluster), function(i){
  #     output[[paste("plot", i, sep="") ]] <- renderPlot({
  #       markergenes <- plotInput()$markergenes
  #       samp.query <- plotInput()$samp.query
  #       ref.cluster <- input$refCluster[i]
  #       ind.ref <- which(colData(newDataRef())$cluster_membership==ref.cluster)
  #       samp.ref <- unname(logcounts(newDataRef())[markergenes, ind.ref])
  #       # FR.test(samp.query[,1:ncol(samp.query)], samp.ref[,1:ncol(samp.ref)],
  #       #         plot.MST=TRUE, label.names = c("Query", "Reference"),
  #       #         main=ref.cluster)
  #       plot(samp.ref[,1], samp.query[,1])
  #     })
  #   })
  # })
  output$temp <- renderPlot({
    ref.cluster <- input$refCluster[1]
    markergenes <- plotInput()$markergenes
    samp.query <- plotInput()$samp.query
          ind.ref <- which(colData(newDataRef())$cluster_membership==ref.cluster)
          samp.ref <- unname(logcounts(newDataRef())[markergenes, ind.ref])
          plot(samp.ref[,1], samp.query[,1])
          # FR.test(samp.query[,1:ncol(samp.query)], samp.ref[,1:ncol(samp.ref)],
          #         plot.MST=TRUE, label.names = c("Query", "Reference"),
          #         main=ref.cluster)
  })
}


# Create Shiny app ----
shinyApp(ui, server)
