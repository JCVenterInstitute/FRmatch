
#' Make data object
#'
#' This function helps to create the data object that can be used in the \code{\link[FRmatch]{FRmatch}} algorithm with input information.
#'
#' @param dat Cell-by-gene expression data in a data frame or equivalent format. First column is cell names, followed by columns for each gene.
#' @param tab Cluster membership of cells in a data frame or equivalent format. First column is cell names, and second column is cluster labels.
#' @param markers A vector of marker genes.
#' @param cluster_marker_info Optionally, a data frame of maker genes for each cluster. See details in \code{\link[FRmatch]{sce.example}}.
#' @param fscores Optionally, a data frame of F-scores for each cluster. See details in \code{\link[FRmatch]{sce.example}}.
#' @param cluster_order Optionally, a vector of ordered cluster names. See details in \code{\link[FRmatch]{sce.example}}.
#'
#' @return A data object of the \link[SingleCellExperiment]{SingleCellExperiment} class
#'
#' @import dplyr
#' @import tibble
#' @importFrom tidyr replace_na
#' @importFrom forcats fct_relevel
#'
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment colData rowData
#' @importMethodsFrom SummarizedExperiment assays assay
# #' @importFrom S4Vectors metadata metadata<-
#'
#' @export
#'

make_data_object <- function(dat, tab, markers,
                             cluster_marker_info=NULL, fscores=NULL, cluster_order=NULL #metadata
){

  ## rename key columns
  names(dat)[1] <- "Sample"
  names(tab) <- c("Sample", "Cluster")

  ## replace special symbols by "_"
  cat("Replace any special symbol in sample and cluster names by '_'. \n")
  dat %<>% mutate(Sample=gsub("-| |\\.|/", "_", Sample))
  tab %<>% mutate(Sample=gsub("-| |\\.|/", "_", Sample), Cluster=gsub("-| |\\.|/", "_", Cluster))

  ## datatable with "Sample", "Cluster", and gene probe columns for constructing sce.object
  dt <- dat %>% inner_join(tab, by="Sample") %>% #inner_join: make sure that cells are in the SAME order!!!
    select("Sample", "Cluster", everything())

  ##----------------------------------##
  ## make SingleCellExperiment object ## sce.object
  ##----------------------------------##

  ## expression matrix
  expr <- dt %>% column_to_rownames("Sample") %>% select(-Cluster) %>% t() #transpose: gene-by-cell
  ## row = genes
  genenames <- rownames(expr) #gene names in order
  ## column = cells
  sampnames <- colnames(expr)

  ## rowData: marker gene
  df_marker_gene <- data.frame("marker_gene"=as.numeric(genenames %in% markers), row.names=genenames)
  ## colData: cluster membership
  df_cluster_membership <- dt %>% column_to_rownames("Sample") %>% select(cluster_membership=Cluster)

  ## make SingleCellExperiment object
  sce.object <- SingleCellExperiment(assays = list(logcounts = expr), #logcounts!!!!!!!!!!!!!
                                     colData = df_cluster_membership,
                                     rowData = df_marker_gene)

  ## metadata:
  metadata(sce.object)$cluster_marker_info <- cluster_marker_info
  metadata(sce.object)$fscores <- fscores
  metadata(sce.object)$cluster_order <- cluster_order

  ## output
  return(sce.object)
}
