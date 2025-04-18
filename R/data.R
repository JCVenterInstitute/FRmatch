#'
#' Example data object
#'
#' This is an example of input data object for \code{\link[FRmatch]{FRmatch}}.
#' The data object is of the \link[SingleCellExperiment]{SingleCellExperiment} class.
#'
#' @format
#' Please follow the naming convention and data formats of this example when constructing your own data object.
#' \describe{
#' \item{\code{assays}}{Data assay \code{logcounts} is the gene-by-cell expression matrix.}
#' \item{\code{rownames}}{Gene names.}
#' \item{\code{rowData}}{Column \code{marker_gene} is a numeric vector indicating if a gene is a marker gene (1) or not (0).}
#' \item{\code{colnames}}{Cell names.}
#' \item{\code{colData}}{Column \code{cluster_membership} is a character vector of cell cluster membership.}
#' \item{\code{metadata}}{List of 3.
#' \itemize{
#'   \item \code{cluster_marker_info} is a data frame of maker genes for each cluster.
#'   \item \code{f_score} is a data frame of F-beta scores for each cluster.
#'   \item \code{cluster_order} is a character vector of ordered cluster names.
#' }}
#' }
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.

"sce.example"
