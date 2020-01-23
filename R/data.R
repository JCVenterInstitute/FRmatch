#'
#' FRmatch input data
#'
#' This page describes the \link[SingleCellExperiment]{SingleCellExperiment} data object that can be used
#' as input data to \code{\link[FRmatch]{FRmatch}}. This example dataset contains 16497 genes and 865 cells.
#'
#' @format
#' Please follow the names and data formats of this example data object when constructing your own data.
#' \describe{
#' \item{\code{metadata}}{List of 3.
#' \itemize{
#'   \item \code{cluster_marker_info} are the cluster-specific NS-Forest marker genes.
#'   \item \code{fsocres} are F-measures associated with the NS-Forest marker genes for each cluster.
#'   \item \code{cluster_order} is some specific order of the clusters to be used for plotting.
#' }}
#' \item{\code{assays}}{Data assay \code{counts} is the expression data matrix (gene by cell).}
#' \item{\code{rownames}}{Gene names.}
#' \item{\code{rowData}}{Column \code{NSF_markers} is a logical vector indicating if a gene is an NS-Forest marker gene.}
#' \item{\code{colnames}}{Cell IDs.}
#' \item{\code{colData}}{Column \code{cluster_membership} is a character vector of cell cluster membership.}
#' }
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} data class.

"sce.example"
