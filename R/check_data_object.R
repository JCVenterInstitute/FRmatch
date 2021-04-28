
#' Check data object
#'
#' We choose to summarize the clustered gene expression data and metadata in a data object of
#' the \code{\link[SingleCellExperiment]{SingleCellExperiment}} class.
#' This function helps to check if the data object is compatible with \code{\link[FRmatch]{FRmatch}}.
#'
#' @param sce.object Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param verbose Boolean variable indicating if to print details to the console. Default: \code{TRUE}.
#'
#' @details For the FR-Match algorithm, the following data elements are essential:
#' \itemize{
#' \item gene expression data matrix
#' \item cell cluster membership
#' \item \emph{informative} marker genes of the reference dataset
#' }
#' In addition, metadata such as \code{cluster_marker_info, f_score, cluster_order} are not essential,
#' but will facilitate visualization and more customized analyses provided in this package.
#'
#' @return A data object that passes this sanity check if no error occurs.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.
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
    stop("'rownames' of this data object is not found. Please see example in help('sce.example').")
  }

  ## rowData
  if(is.null(rowData(sce.object)$marker_gene)){
    stop("'marker_gene' is not found in the rowData of this data object. Please see example in help('sce.example').")
  }

  #############
  ## colData ##
  #############
  ## colnames
  if(is.null(colnames(sce.object))){
    stop("'colnames' of this data object is not found. Please see example in help('sce.example').")
  }
  ## colData
  if(is.null(colData(sce.object)$cluster_membership)){
    stop("'cluster_membership' is not found in the colData of this data object. Please see example in help('sce.example').")
  }

  ##############
  ## metadata ##
  ##############
  ## cluster_marker_info
  if(is.null(sce.object@metadata$cluster_marker_info)){
    if(verbose) cat("'cluster_marker_info' is not available. \n")
  }
  ## f_score
  if(is.null(sce.object@metadata$f_score)){
    if(verbose) cat("'f_score' is not available. \n")
  }
  ## cluster_order
  if(is.null(sce.object@metadata$cluster_order)){
    sce.object@metadata$cluster_order <- sort(unique(colData(sce.object)$cluster_membership))
    if(verbose) cat("'cluster_order' (alphabetical order) is added to the metadata of this data object. \n")
  }

  ######################
  ## final sce.object ##
  ######################
  if(verbose) cat("Data object check complete. \n")
  return(sce.object)
}








