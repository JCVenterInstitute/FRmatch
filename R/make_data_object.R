
#' Make data object
#'
#' This function creates the data object used in \code{\link[FRmatch]{FRmatch}}.
#' See an example of the data object in \code{\link[FRmatch]{sce.example}}.
#'
#' @param dat Cell-by-gene expression data in a data frame or equivalent format. The first column is cell names, followed by columns of each gene.
#' @param tab Cluster membership of each cell in a data frame or equivalent format. The first column is cell names, and the second column is cluster labels.
#' @param markers A vector of marker genes.
#' @param cluster_marker_info Optionally, a data frame of maker genes for each cluster. \code{marker} is the unique set of marker genes for all clusters.
#' @param f_score Optionally, a data frame of F-beta scores for each cluster.
#' @param cluster_order Optionally, a vector of ordered cluster names, which may reflect the taxonomy of the cell type clusters.
#'
#' @return A data object of the \link[SingleCellExperiment]{SingleCellExperiment} class
#'
#' @import dplyr
#' @import tibble
#' @importFrom magrittr %<>%
#' @importFrom tidyr replace_na
#' @importFrom forcats fct_relevel
#'
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment colData rowData
#' @importMethodsFrom SummarizedExperiment assays assay
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export
#'

make_data_object <- function(dat, tab, markers,
                             cluster_marker_info=NULL, f_score=NULL, cluster_order=NULL #metadata
){

  ## rename key columns
  names(dat)[1] <- "Sample"
  names(tab) <- c("Sample", "Cluster")

  ## replace special symbols by "_"
  # cat("Replace any special symbol in sample and cluster names by '_'. \n")
  # dat %<>% mutate(Sample=gsub("-| |\\.|/", "_", Sample))
  # names(dat) <- gsub("-| |\\.|/", "_", names(dat))
  # tab %<>% mutate(Sample=gsub("-| |\\.|/", "_", Sample), Cluster=gsub("-| |\\.|/", "_", Cluster))
  # markers <- gsub("-| |\\.|/", "_", markers)
  # if(!is.null(cluster_marker_info)){
  #   cluster_marker_info %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName), markerGene=gsub("-| |\\.|/", "_", markerGene))
  # }
  # if(!is.null(f_score)){
  #   f_score %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName))
  # }
  # if(!is.null(cluster_order)){
  #   cluster_order <- gsub("-| |\\.|/", "_", cluster_order)
  # }

  ## data table with "Sample", "Cluster", and gene columns for constructing the sce.object
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
  if(!is.null(cluster_marker_info)){
    names(cluster_marker_info) <- c("clusterName", "markerGene")
    # cluster_marker_info %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName))
    # cluster_marker_info %<>% mutate(markerGene=gsub("-| |\\.|/", "_", markerGene))
  }
  if(!is.null(f_score)){
    names(f_score) <- c("clusterName", "score")
    # f_score %<>% mutate(clusterName=gsub("-| |\\.|/", "_", clusterName))
  }
  if(!is.null(cluster_order)){
    # cluster_order <- gsub("-| |\\.|/", "_", cluster_order)
  }
  metadata(sce.object)$cluster_marker_info <- cluster_marker_info
  metadata(sce.object)$f_score <- f_score
  metadata(sce.object)$cluster_order <- cluster_order

  ## output
  return(sce.object)
}
