library(shiny)
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

      ## Input: Specify random seed ----
      numericInput("seed", "Please set a random seed:", 999),

      ## Input: Specify spliting fraction ----
      sliderInput("splitFrac", "Data subsampling fraction:",
                  0.2, 0.8, 0.5, step=0.1),
      ## Include clarifying text ----
      helpText("Note: the above fraction of cells will be selected",
               "in proportion to the cluster sizes from the reference dataset,",
               "and the rest of the cells will be selected from the query dataset."),

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
                 plotOutput("dropouts", height="600px")),
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
                 helpText("Please ignore the error message and",
                          "choose one or more reference cluster(s)."),
                 uiOutput("plotMST"))

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
  newDataRef <- eventReactive(input$updateButton, {
    newData()$newsce.ref
  }, ignoreNULL = FALSE)
  newDataQuery <- eventReactive(input$updateButton, {
    newData()$newsce.query
  }, ignoreNULL = FALSE)

  ## run FRmatch
  results <- eventReactive(input$updateButton, {
    FRmatch(newDataQuery(), newDataRef())
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
    check_dropout(newDataRef(), return.value=TRUE, plot.dropout=TRUE)
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
  output$MST <- renderPlot({
    markergenes <- newDataRef()@metadata$cluster_marker_info$markerGene
    query.cluster <- input$queryCluster
    ind.query <- which(colData(newDataQuery())$cluster_membership==query.cluster)
    samp.query <- unname(logcounts(newDataQuery())[markergenes, ind.query])
    par(mfrow=c(ceiling(length(input$refCluster)/2),2))
    for(ref.cluster in input$refCluster){
      ind.ref <- which(colData(newDataRef())$cluster_membership==ref.cluster)
      samp.ref <- unname(logcounts(newDataRef())[markergenes, ind.ref])
      # samp.ref
      FR.test(samp.query[,1:ncol(samp.query)], samp.ref[,1:ncol(samp.ref)],
              plot.MST=TRUE, label.names = c("Query", "Reference"),
              main=ref.cluster)
    }
  })
  # plotHeight <- reactive(350 * ceiling(length(input$refCluster)/2))
  output$plotMST <- renderUI({
    plotOutput("MST", height = 300*ceiling(length(input$refCluster)/2))
  })

}


# Create Shiny app ----
shinyApp(ui, server)
