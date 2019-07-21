
#' Sanity check for input data object
#'
#' For input data, we choose to summarize the gene expression data and metadata needed for \code{FRmatch} in
#' the \code{\link[SingleCellExperiment]{SingleCellExperiment}} class, which is a light-weighted container for single-cell genomics data.
#' This function helps to check if users' data objects fit the \code{FRmatch} context.
#'
#' @param sce.object A \code{SingleCellExperiment} object filled with necessary information for \code{FRmatch}.
#' See details in \code{\link[FRmatch]{sce.example}}.
#'
#' @details For FRmatch, the following data items are essential:
#' \itemize{
#' \item a gene expression data matrix (gene by cell)
#' \item cluster membership information for each cell
#' \item NS-Forest (or your own) marker genes for the reference clusters
#' }
#' Additional items, such as F-meansure (i.e. a quality score associated with NS-Forest markers) and preferred cluster order,
#' are not essential, but will facilitate graphical tools provided in this package.
#'
#' @return A \code{SingleCellExperiment} object that passes this sanity check if no error occurs.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.
#' @export


check_data <- function(sce.object){

  ###############
  ## metatdata ##
  ###############
  ## cluster_marker_info
  if(is.null(sce.object@metadata$cluster_marker_info)){
    stop("@metadata$cluster_marker_info is not found. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }
  if(!is.null(sce.object@metadata$cluster_marker_info)){
    if(is.null(sce.object@metadata$cluster_marker_info$cluster)){
      stop("@metadata$cluster_marker_info$cluster is not found. Please see example in help('sce.example'). \n
           Data object check not complete. \n")
    }
    if(is.null(sce.object@metadata$cluster_marker_info$markerGene)){
      stop("@metadata$cluster_marker_info$markerGene is not found. Please see example in help('sce.example'). \n
           Data object check not complete. \n")
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
  if(!"logcounts" %in% names(assays(sce.object))){
    stop("'logcounts' not in names(assays(<SingleCellExperiment>)). Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }

  #############
  ## rowData ##
  #############
  ## rownames
  if(is.null(rownames(sce.object))){
    stop("rownames of this data object is not found. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
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
    stop("colnames of this data object is not found. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }
  ## colData
  if(is.null(SingleCellExperiment::colData(sce.object)$cluster_membership)){
    stop("'cluster_membership' is not found in colData. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }

  ## check cluster names
  if(length(base::setdiff(unique(SingleCellExperiment::colData(sce.object)$cluster_membership), sce.object@metadata$cluster_order))>0 |
     length(base::setdiff(sce.object@metadata$cluster_order, unique(SingleCellExperiment::colData(sce.object)$cluster_membership)))>0){
    stop("Cluster names in colData and metadata do not match. Please use consistent cluster names. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
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
  if(!"logcounts" %in% names(assays(sce.object))){
    stop("'logcounts' not in names(assays(<SingleCellExperiment>)). Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }

  #############
  ## rowData ##
  #############
  ## rownames
  if(is.null(rownames(sce.object))){
    stop("rownames of this data object is not found. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }

  #############
  ## colData ##
  #############
  ## colnames
  if(is.null(rownames(sce.object))){
    stop("colnames of this data object is not found. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }
  ## colData
  if(is.null(SingleCellExperiment::colData(sce.object)$cluster_membership)){
    stop("'cluster_membership' is not found in colData. Please see example in help('sce.example'). \n
         Data object check not complete. \n")
  }

  ###############
  ## metatdata ##
  ###############

  ## cluster_order
  if(is.null(sce.object@metadata$cluster_order)){
    sce.object@metadata$cluster_order <- sort(unique(SingleCellExperiment::colData(sce.object)$cluster_membership))
    cat("@metadata$cluster_order was not available originally. It is replaced by alphabetical order of unique cluster names. \n")
  }

  ######################
  ## final sce.object ##
  ######################
  cat("Data object check complete. \n")
  return(sce.object)
  }






