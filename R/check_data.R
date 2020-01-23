
#' Sanity check of input data object
#'
#' For input data, we choose to summarize the gene expression data and metadata needed for \code{\link[FRmatch]{FRmatch}} in a data object of
#' the \code{\link[SingleCellExperiment]{SingleCellExperiment}} class, which is a light-weighted container for single-cell genomics data.
#' This function helps to check if users' data objects can be used in the \code{FRmatch} function.
#'
#' @param sce.object A \code{FRmatch} input data object. See details in \code{\link[FRmatch]{sce.example}}.
#'
#' @details For the FR-Match method, the following data are essential:
#' \itemize{
#' \item a gene expression data matrix (gene by cell)
#' \item cell cluster membership information
#' \item NS-Forest (or your own) marker genes for the reference dataset
#' }
#' In addition, information such as F-meansure and cluster order are not essential, but will facilitate visualization and customized analysis
#' provided in this package.
#'
#' @return A data object that passes this sanity check if no error occurs.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} data class.
#'
#' @import dplyr
#' @import methods
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importMethodsFrom SingleCellExperiment colData rowData
#' @importMethodsFrom SummarizedExperiment assays
#'
#' @export


check_data <- function(sce.object){

  ##############
  ## metadata ##
  ##############
  ## cluster_marker_info
  if(is.null(sce.object@metadata$cluster_marker_info)){
    stop("@metadata$cluster_marker_info is not found. Please see example in help('sce.example').")
  }
  if(!is.null(sce.object@metadata$cluster_marker_info)){
    if(is.null(sce.object@metadata$cluster_marker_info$cluster)){
      stop("@metadata$cluster_marker_info$cluster is not found. Please see example in help('sce.example').")
    }
    if(is.null(sce.object@metadata$cluster_marker_info$markerGene)){
      stop("@metadata$cluster_marker_info$markerGene is not found. Please see example in help('sce.example').")
    }
  }
  ## fscores
  if(is.null(sce.object@metadata$fscores)){
    cat("@metadata$fscores is not available. \n")
  }
  ## cluster_order
  if(is.null(sce.object@metadata$cluster_order)){
    sce.object@metadata$cluster_order <- sort(unique(sce.object@metadata$cluster_marker_info$cluster))
    cat("@metadata$cluster_order was not available originally. It is replaced by alphabetical order of unique cluster names. \n")
  }

  ############
  ## assays ##
  ############
  if(!"counts" %in% names(assays(sce.object))){
    stop("'counts' not in names(assays(sce.object)). Please see example in help('sce.example').")
  }

  #############
  ## rowData ##
  #############
  ## rownames
  if(is.null(rownames(sce.object))){
    stop("rownames of this data object is not found. Please see example in help('sce.example').")
  }
  ## check gene names
  if(length(base::setdiff(sce.object@metadata$cluster_marker_info$markerGene,rownames(sce.object)))>0){
    warning("At least 1 marker genes not presented in the rownames of this data object. Please make sure consistent names of genes are used. Pay attention to special symbols such as '-', '.', '_'. \n")
  }
  ## rowData
  if(is.null(rowData(sce.object)$NSF_markers)){
    rowData(sce.object)$NSF_markers <- rownames(sce.object) %in% unique(sce.object@metadata$cluster_marker_info$markerGene)
    cat("'NSF_markers' column is added to rowData. \n")
  }
  ## for updating after filter.cluster??
  # if(!identical(rownames(sce.object)[rowData(sce.object)$NSF_markers], sce.object@metadata$cluster_marker_info$markerGene)){
  #   rowData(sce.object)$NSF_markers <- rownames(sce.object) %in% unique(sce.object@metadata$cluster_marker_info$markerGene)
  #   cat("'NSF_markers' column is updated according to metadata. \n")
  # }

  #############
  ## colData ##
  #############
  ## colnames
  if(is.null(rownames(sce.object))){
    stop("colnames of this data object is not found. Please see example in help('sce.example').")
  }
  ## colData
  if(is.null(colData(sce.object)$cluster_membership)){
    stop("'cluster_membership' is not found in colData. Please see example in help('sce.example').")
  }

  ## check cluster names
  if(length(base::setdiff(unique(colData(sce.object)$cluster_membership), sce.object@metadata$cluster_order))>0 |
     length(base::setdiff(sce.object@metadata$cluster_order, unique(colData(sce.object)$cluster_membership)))>0){
    stop("Cluster names in colData and metadata do not match. Please use consistent cluster names. Please see example in help('sce.example').")
  }

  ######################
  ## final sce.object ##
  ######################
  cat("Data object check complete. \n")
  return(sce.object)
}

#######################################################################################################


check_query_data <- function(sce.object){

  ############
  ## assays ##
  ############
  if(!"counts" %in% names(assays(sce.object))){
    stop("'counts' not in names(assays(sce.object)). Please see example in help('sce.example').")
  }

  #############
  ## rowData ##
  #############
  ## rownames
  if(is.null(rownames(sce.object))){
    stop("rownames of this data object is not found. Please see example in help('sce.example').")
  }

  #############
  ## colData ##
  #############
  ## colnames
  if(is.null(rownames(sce.object))){
    stop("colnames of this data object is not found. Please see example in help('sce.example').")
  }
  ## colData
  if(is.null(colData(sce.object)$cluster_membership)){
    stop("'cluster_membership' is not found in colData. Please see example in help('sce.example').")
  }

  ###############
  ## metatdata ##
  ###############

  ## cluster_order
  if(is.null(sce.object@metadata$cluster_order)){
    sce.object@metadata$cluster_order <- sort(unique(colData(sce.object)$cluster_membership))
    cat("@metadata$cluster_order was not available originally. It is replaced by alphabetical order of unique cluster names. \n")
  }

  ######################
  ## final sce.object ##
  ######################
  cat("Data object check complete. \n")
  return(sce.object)
  }






