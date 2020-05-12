
#' Check data object
#'
#' We choose to summarize the clustered gene expression data and metadata in a data object of
#' the \code{\link[SingleCellExperiment]{SingleCellExperiment}} class.
#' This function helps to check if the data object is compatible with \code{\link[FRmath]{FRmatch}}.
#'
#' @param sce.object Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param verbose Boolean variable indicating if to print details to the console. Default: \code{TRUE}.
#'
#' @details For the FR-Match algorithm, the following data elements are essential:
#' \itemize{
#' \item a gene expression data matrix
#' \item cell cluster membership
#' \item \emph{informative} marker genes of the reference dataset
#' }
#' In addition, metadata such as \code{cluster_marker_info, fscores, cluster_order} are not essential,
#' but will facilitate visualization and more customized analyses provided in this package.
#'
#' @return A data object that passes this sanity check if no error occurs.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} data class.
#'
#' @import dplyr
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment colData rowData
#' @importMethodsFrom SummarizedExperiment assays
#' @importFrom S4Vectors metadata metadata<-
#'
#' @export


check_data_object <- function(sce.object, verbose=TRUE){

  ############
  ## assays ##
  ############
  if(!"logcounts" %in% names(assays(sce.object))){
    stop("'logcounts' is not found in the data object. Please see example in help('sce.example').")
  }

  #############
  ## rowData ##
  #############
  ## rownames
  if(is.null(rownames(sce.object))){
    stop("rownames of this data object is not found. Please see example in help('sce.example').")
  } else {
    if(verbose) cat("Replace any special character in rownames by '_'. \n")
    rownames(sce.object) <- gsub("-| |\\.|/", "_", rownames(sce.object))
  }
  ## rowData
  if(is.null(rowData(sce.object)$marker_gene)){
    stop("'marker_gene' is not found in the rowData of this data object. Please see example in help('sce.example').")
    # rowData(sce.object)$marker_gene <- as.numeric(rownames(sce.object) %in% unique(sce.object@metadata$cluster_marker_info$markerGene))
    # cat("'marker_gene' column is added to rowData of this data object. \n")
  }

  # ## check gene names
  # if(length(base::setdiff(sce.object@metadata$cluster_marker_info$markerGene,rownames(sce.object)))>0){
  #   warning("At least 1 marker genes not presented in the rownames of this data object. Please make sure consistent names of genes are used. Pay attention to special symbols such as '-', '.', '_'. \n")
  # }
  # ## for updating after filter.cluster??
  # if(!identical(rownames(sce.object)[rowData(sce.object)$NSF_markers], sce.object@metadata$cluster_marker_info$markerGene)){
  #   rowData(sce.object)$NSF_markers <- rownames(sce.object) %in% unique(sce.object@metadata$cluster_marker_info$markerGene)
  #   cat("'NSF_markers' column is updated according to metadata. \n")
  # }

  #############
  ## colData ##
  #############
  ## colnames
  if(is.null(colnames(sce.object))){
    stop("colnames of this data object is not found. Please see example in help('sce.example').")
  } else {
    if(verbose) cat("Replace any special character in colnames by '_'. \n")
    colnames(sce.object) <- gsub("-| |\\.|/", "_", colnames(sce.object))
  }
  ## colData
  if(is.null(colData(sce.object)$cluster_membership)){
    stop("'cluster_membership' is not found in the colData of this data object. Please see example in help('sce.example').")
  }

  # ## check cluster names
  # if(length(base::setdiff(unique(colData(sce.object)$cluster_membership), sce.object@metadata$cluster_order))>0 |
  #    length(base::setdiff(sce.object@metadata$cluster_order, unique(colData(sce.object)$cluster_membership)))>0){
  #   stop("Cluster names in colData and metadata do not match. Please use consistent cluster names. Please see example in help('sce.example').")
  # }

  ##############
  ## metadata ##
  ##############
  ## cluster_marker_info
  if(is.null(sce.object@metadata$cluster_marker_info)){
    if(verbose) cat("'cluster_marker_info' is not available. \n")
  }
  ## fscores
  if(is.null(sce.object@metadata$fscores)){
    if(verbose) cat("'fscores' is not available. \n")
  }
  ## cluster_order
  if(is.null(sce.object@metadata$cluster_order)){
    sce.object@metadata$cluster_order <- sort(unique(sce.object@metadata$cluster_marker_info$cluster))
    if(verbose) cat("'cluster_order' (alphabetical order) is added to the metadata of this data object. \n")
  }

  ######################
  ## final sce.object ##
  ######################
  if(verbose) cat("Data object check complete. \n")
  return(sce.object)
}








