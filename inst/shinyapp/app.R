library(shiny)
library(shinycssloaders)
library(shinythemes)

# library(plotly)
library(dplyr)
library(tibble)
library(FRmatch)
library(SingleCellExperiment)


load("../data/sce-layer1-15clusters.rda")
load("../data/sce-layer1-topNodes.rda")

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
## UI
#######################################################################################################

# Define UI for dataset viewer app ----
ui <- fluidPage(

  # Shiny theme
  theme = shinytheme("cosmo"),

  # App title ----
  titlePanel(div("FR-Match Demo",
                 img(height = 35, width = 200,
                     src = "JCVI-Logo-Inline-Black.png",
                     class = "pull-right"))),

  # Sidebar layout with input and output definitions ----
  sidebarLayout(

    # Sidebar panel for inputs ----
    sidebarPanel(

      h3("Pre-loaded data"),
      p("In this Shiny App, we preloaded two example datasets to start with. These two datasets are essentiall
        the same data, but with different cell type cluster labels. The 'Layer1' data is from Layer 1 of the
        human middle temporal gyrus brain region, with 15 cell types defined in Boldog et al. (2018). The
        'Layer1-topNode' groups those cell types to the top level broad class of brain cells."),

      ## Input: Select datasets ----
      selectInput("querydata", "Choose a query dataset:",
                  choices = c("Layer1", "Layer1-topNodes")),
      selectInput("refdata", "Choose a reference dataset:",
                  choices = c("Layer1-topNodes", "Layer1")),
      helpText("Boldog, Eszter, et al. 'Transcriptomic and morphophysiological evidence for a specialized human
               cortical GABAergic cell type.' Nature neuroscience 21.9 (2018): 1185-1195."),

      # ## Input: Checkbox if to impute reference ----
      # checkboxInput("imputation", "Imputation", FALSE),
      # ## Include clarifying text ----
      # helpText(p("Note: by checking the Imputation box,",
      #            "FRmatch will impute the dropout values for each marker gene",
      #            em("only"), "in the cluster that it marks and",
      #            em("only"), "in the reference dataset.")),

      # Horizontal line ----
      tags$hr(),

      h3("Cross-validation data"),
      p("For illustration purpose, we demonstrate a cross-validation on the pre-loaded data. Use the following
        to split the cells into the query and reference datasets."),

      ## Input: Specify random seed ----
      numericInput("seed", "Please set a random seed:", 100),

      ## Input: Specify spliting fraction ----
      sliderInput("splitFrac", "Data splitting fraction:",
                  0.2, 0.8, 0.5, step=0.1),
      ## Include clarifying text ----
      helpText("The above fraction of cells will be selected in proportion to the cluster sizes
               from the reference dataset, and the rest of the cells will be selected from the query dataset."),

      # Horizontal line ----
      tags$hr(),

      ## Input: actionButton() to defer the rendering of output ----
      actionButton("updateButton", "Run FRmatch",
                   class = "btn-primary")

    ), #close sidebarPanel


    # Main panel for displaying outputs ----
    mainPanel(

      tabsetPanel(
        tabPanel("Data",
                 h2("Input data"),
                 p("A quick view of the data class, data dimensions, row data, column data, metatdata, etc."),
                 h3("Query dataset"),
                 verbatimTextOutput("viewQuery"),
                 tableOutput("tableQuery"),
                 h3("Reference dataset"),
                 verbatimTextOutput("viewRef"),
                 tableOutput("tableRef")),
        # tabPanel("Dropouts",
        #          h2("Check %expressed in the reference data"),
        #          checkboxInput("afterImputation", "Show after imputation", FALSE),
        #          plotOutput("dropouts", height="auto")),
        tabPanel("Cluster size",
                 h2("Overview of clusters"),
                 p("FR-Match is a cluster-level matching algorithm. It has many tuning parameters that
                   are associated with cluster sizes, such as the filtering step of small clusters, and
                   the subsampling size uder the iterative procedure. Here, we provide an overview of the
                   clusters and a simple comparison of their sizes in the query and reference datasets."),
                 plotOutput("clusterSize")),
        tabPanel("Bracode plot",
                 h2("'Barcoding' clusters by marker genes"),
                 p("FR-Match uses informative marker genes as a dimensionality reduction tool that is key
                   to the matching performance. We may utilize the `barcode` plot to get insights of how well
                   are the marker genes tagging the reference clusters. Intuitively, we may think of these plots
                   as the scanning barcode of products in a grocery store."),
                 uiOutput("clusterSelection"),
                 uiOutput("plotBarcode")),
        tabPanel("MST plot",
                 h2("Minimum Spanning Tree"),
                 p("The core of FR-Match is a graphical model based on the Minimum Spanning Tree (MST). An great
                   advantage of using such a graphial model is that we may visualize the data cloud and visually
                   exam the relationships of the cells in query and referenct clusters on the MST plot.
                   An interwoven MST suggests a match of the query and reference clusters."),
                 uiOutput("querySelection"),
                 uiOutput("refSelection"),
                 helpText("Please ignore the error message and",
                          "choose one or more reference cluster(s)."),
                 uiOutput("plotMST")),
        tabPanel("FR-Match results",
                 h3("Recommended matches"),
                 plotOutput("matches") %>% withSpinner(color="#0dc5c1"),
                 sliderInput("sigLvl", "Significance level:",
                             0, 0.2, 0.05, step=0.01, width="400px"),
                 helpText("By setting smaller significance level, there will be more matches found."),
                 hr(),
                 h3("Distribution of adjusted p-values"),
                 plotOutput("padj") %>% withSpinner(color="#0dc5c1"))

      ) #close tabsetPanel

    ) #closde mainPanel

  ) #close sidebarLayout

) #close fliudPage

#######################################################################################################
## server
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

  ## splitting data
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
    FRmatch(newDataQuery(), newDataRef(), imputation=input$imputation, filter.size=0, subsamp.iter=101)
  }, ignoreNULL = FALSE)

  ##------ tab: Data ------##

  ## look at query data
  output$viewQuery <- renderPrint({
    newDataQuery()
  })
  # output$tableQuery <- renderTable({
  #   tab.query <- table(colData(queryInput())$cluster_membership)
  #   tab.newquery <- table(colData(newDataQuery())$cluster_membership)
  #   table.query <- cbind(tab.query, tab.newquery) %>% data.frame() %>% rownames_to_column()
  #   colnames(table.query) <- c("Cluster", "Size", "Selected")
  #   table.query
  # }, rownames=TRUE)

  ## look at ref data
  output$viewRef <- renderPrint({
    newDataRef()
  })
  # output$tableRef <- renderTable({
  #   tab.ref <- table(colData(refInput())$cluster_membership)
  #   tab.newref <- table(colData(newDataRef())$cluster_membership)
  #   table.ref <- cbind(tab.ref, tab.newref) %>% data.frame() %>% rownames_to_column()
  #   colnames(table.ref) <- c("Cluster", "Size", "Selected")
  #   table.ref
  # }, rownames=TRUE)

  ##------ tab: Dropouts ------##

  # ## dropout plot
  # output$dropouts <- renderPlot({
  #   if(input$afterImputation){
  #     newDataRefImputation <- FRmatch:::impute.zero(newDataRef())
  #     plot_nonzero(newDataRefImputation, return.value=FALSE, return.plot=TRUE)
  #   }
  #   else plot_nonzero(newDataRef(), return.value=FALSE, return.plot=TRUE)
  # }, height = function() {
  #   session$clientData$output_dropouts_width
  # })

  ##------ tab: Cluster size ------##

  ## cluster size plot
  output$clusterSize <- renderPlot({
    plot_clusterSize(newDataQuery(), newDataRef(), name.E1 = "Layer1", name.E2 = "Layer1-topNodes")
  }, height = 800)

  ##------ tab: Barcode plot ------##

  ## barcode plot
  output$clusterSelection <- renderUI({
    selectInput("cluster", "Select reference cluster:", choices = newDataRef()@metadata$cluster_order,
                multiple = FALSE)
  })
  output$barcode <- renderPlot({
    plot_cluster_by_markers(newDataRef(), cluster.name = input$cluster, name.E1 = "ref")
    # par(mfrow=c(ceiling(length(input$cluster)/2),2))
    # for(cluster in input$cluster){
    #   plot_cluster_by_markers(newDataRef(), cluster.name = cluster, name.E1 = "ref")
    # }
  })
  output$plotBarcode <- renderUI({
    plotOutput("barcode")
    # plotOutput("barcode", height = 500*ceiling(length(input$cluster)/2))
  })

  ##------ tab: FR-Match results ------##

  ## matches plot
  output$matches <- renderPlot({
    plot_FRmatch(results(), sig.level=input$sigLvl, reorder=FALSE)
  })

  ## p-values plot
  output$padj <- renderPlot({
    plot_FRmatch(results(), type="padj", sig.level=input$sigLvl, reorder=FALSE)
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

#######################################################################################################
#######################################################################################################

# Create Shiny app ----
shinyApp(ui, server)

