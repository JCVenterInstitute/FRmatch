#'
#' Details about input data
#'
#' This page documents the \code{SingleCellExperiment} object and necessary information needed for \code{\link[FRmatch]{FRmatch}}
#' and other functions in this package.
#' DESCRIPTION OF RAW DATA. This example scRNA-seq data contain 865 cells with 16497 genes.
#'
#' @format
#' Please follow the naming convention and format of this example.
#' \describe{
#' \item{metadata}{List of 3.
#' \itemize{
#'   \item \code{cluster_marker_info} contains NS-Forest marker genes of each clusters.
#'   \item \code{fsocres} is the F-measure classification score for the set of NS-Forest markers for each cluster.
#'   \item \code{cluster_order} is some specific order of clusters to facilitate plot.
#' }}
#' \item{assays}{One data assay \code{logcounts}, which is an expression data matrix (gene by cell) with log-transformed counts.}
#' \item{rownames}{Gene names.}
#' \item{rowData}{One column \code{NSF_markers}, which is a logical vector indicating if a gene is an NS-Forest marker gene.}
#' \item{colnames}{Cell IDs.}
#' \item{colData}{One column \code{cluster_membership}, which is a character vector of cluster membership of each cell.}
#' }
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.

"sce.example"
