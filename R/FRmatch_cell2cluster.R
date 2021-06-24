
#' \code{FRmatch} cell-to-cluster extension
#'
#' This function is an extension of the original \code{\link[FRmatch]{FRmatch}} to assign each query cell with a reference cluster label.
#' Please see Details for the extension.
#'
#' @param sce.query Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for query experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param sce.ref Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for reference experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
# #' @param imputation INACTIVE. Logical variable indicating if to impute expression zero values for the reference experiment. Default: \code{FALSE}.
# #' See details in \code{\link[FRmatch]{impute_dropout}}.
#' @param filter.size,filter.fscore Filtering small/poor-quality clusters. Default: \code{filter.size=10}, filter based on the number
#' of cells per cluster; \code{filter.fscore=NULL}, filter based on the F-beta score associated with the cell cluster if available (numeric value).
#' @param subsamp.size,subsamp.iter,subsamp.seed Cluster size, number of iterations, and random seed.
#' Default: \code{10, 1001, 1}, respectively.
#' @param numCores Number of cores for parallel computing.
#' Default: \code{NULL}, use the maximum number of cores detected by \code{\link[parallel]{detectCores}} if not specified (an integer).
#' @param prefix Prefix names for query and reference clusters. Default: \code{prefix=c("query.", "ref.")}.
#' @param verbose Numeric value indicating levels of details to be printed. Default: \code{1}, only print major steps.
#' If \code{0}, no verbose; if \code{2}, print all.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FRtest}}.
#'
#' @return A data frame with columns:
#' \item{cell}{Query cell ID.}
#' \item{cluster}{Cluster membership of query cells.}
#' \item{match.cell2cluster}{Matched reference cluster for the query cell by \code{FRmatch_cell2cluster}.}
#' \item{rmax.cell2cluster}{Row maximum of \code{pmat.cell2cluster} (see below).}
#' And concatenated by the columns from the matrix:
#' \item{pmat.cell2cluster}{Cell-by-cluster (a.k.a. query cell by reference cluster) matrix of p-values by \code{FRmatch_cell2cluster}.}
#'
#' @details
#' Apply \code{FRmatch} with its iterative subsampling scheme, which is a bootstrap-like approach to randomly select a smaller set of cells
#' from a query cluster and quantify the confidence score of the selected cells belonging to certain reference cell type
#' using the p-value outputted from \code{FRmatch} (i.e. a larger p-value indicates higher probability of a match, and vice versa).
#'
#' Assign the cluster-level p-value to each selected query cell, and updated the assigned p-value if the query cell is reselected
#' from the iterative procedure and assigned a higher p-value. The output from the cell-to-cluster extension is a cell-by-cluster
#' (a.k.a. query cell by reference cluster) matrix of p-values.
#'
#' @author Yun Zhang, \email{zhangy@jcvi.org};
#' Brian Aevermann, \email{baeverma@jcvi.org};
#' Richard Scheuermann, \email{RScheuermann@jcvi.org}.
#'
# @seealso Visualization of matching results using \code{\link[FRmatch]{plot_FRmatch}}, \code{\link[FRmatch]{plot_bi_FRmatch}}.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' FRmatch_cell2cluster(sce.example, sce.example)
#' }
#'
#' @export

FRmatch_cell2cluster <- function(sce.query, sce.ref, #imputation=FALSE,
                                 filter.size=10, filter.fscore=NULL, #filtering clusters
                                 subsamp.size=10, subsamp.iter=1001, subsamp.seed=1, #subsampling
                                 numCores=NULL, prefix=c("query.", "ref."),
                                 verbose=1, ...){

  ## check data object
  if(verbose>0) cat("* Check query data object. \n")
  sce.query <- check_data_object(sce.query, verbose=(verbose>1))
  if(verbose>0) cat("* Check reference data object. \n")
  sce.ref <- check_data_object(sce.ref, verbose=(verbose>1))

  ## filtering small or low fscore clusters
  if(verbose>0) cat("* Filtering small clusters: query and reference clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)){
    if(verbose>0) cat("* Filtering low F-beta score clusters: reference cluster with F-beta score <", filter.fscore, " are not considered. \n")
  }
  sce.query <- filter_cluster(sce.query, filter.size=filter.size, filter.fscore=NULL)
  sce.ref <- filter_cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)

  ## imputation
  # if(imputation){
  #   sce.ref <- impute_zero(sce.ref)
  #   cat("Imputation is applied. \n")
  # }

  ## extract info from sce.objects
  querydat <- assay(sce.query) #matrix
  refdat <- assay(sce.ref)
  membership.query <- colData(sce.query)$cluster_membership
  membership.ref <- colData(sce.ref)$cluster_membership
  order.query <- sce.query@metadata$cluster_order
  order.ref <- sce.ref@metadata$cluster_order

  ## dimension reduction by selecting only marker genes
  markergenes <- rownames(sce.ref)[rowData(sce.ref)$marker_gene==1]
  markergenes.common <- intersect(markergenes, rownames(sce.query))
  querydat.reduced <- querydat[markergenes.common,]
  refdat.reduced <- refdat[markergenes.common,]
  if(verbose>0) cat("* Feature selection:", length(markergenes.common), "out of", length(markergenes),
                    "reference marker genes are presented in the query experiment. \n")
  if(length(markergenes.common)!=length(markergenes)){
    if(verbose>1) cat(paste(setdiff(markergenes,markergenes.common),collapse = ", "),"not found. \n")
    # warning(paste(setdiff(markergenes,markergenes.common),collapse = ", ")," not found.")
  }

  ## split data by cluster membership and store in a list
  datlst.query <- lapply(split(as.data.frame(t(querydat.reduced)), membership.query),t)
  datlst.query <- datlst.query[order.query]
  clusterNames.query <- paste0(prefix[1],names(datlst.query))
  names(datlst.query) <- clusterNames.query
  datlst.ref <- lapply(split(as.data.frame(t(refdat.reduced)), membership.ref),t)
  datlst.ref <- datlst.ref[order.ref]
  clusterNames.ref <- paste0(prefix[2],names(datlst.ref))
  names(datlst.ref) <- clusterNames.ref

  ncluster.query <- length(datlst.query)
  ncluster.ref <- length(datlst.ref)
  if(verbose>0) cat("* Comparing", ncluster.query, "query clusters with", ncluster.ref, "reference clusters... \n")

  ###########################################
  ## FR comparison between two experiments ##
  ###########################################

  ##------ USE PARALLEL COMPUTING!!! ------##
  ## use all available cores if not specified
  if(is.null(numCores)) numCores <- detectCores()
  if(verbose>0) cat("** Parallel computing with", numCores, "cores. \n")

  ## prepare data for each combination of query cluster and ref cluster pair
  paired.datlst.query <- rep(datlst.query, each=ncluster.ref)
  paired.datlst.ref <- rep(datlst.ref, times=ncluster.query)
  ##---------------------------------------##

  if(verbose>0) cat("** method = cell2cluster", "| subsamp.size =", subsamp.size, "| subsamp.iter =", subsamp.iter, "\n")

  results <- mcmapply(
    function(samp1,samp2){
      set.seed(subsamp.seed)
      FRtest_cell2cluster(samp1, samp2, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, ...)
    },
    paired.datlst.query, paired.datlst.ref,
    mc.cores = numCores)

  if(verbose>0) cat("* Done parallel and returning results. \n")

  ## reformatting results from mcmapply
  names(results) <- rep(clusterNames.ref, times=ncluster.query) #for colnames in pmat.cell2cluster
  f <- rep(names(datlst.query), each=ncluster.ref)
  rstlst.query <- split(results, f)[clusterNames.query] #make sure the list is ordered by query cluster order
  rstlst.query.pmat <- lapply(rstlst.query, function(z) do.call("cbind",z))
  pmat.cell2cluster <- do.call("rbind", rstlst.query.pmat)
  if(sum(is.na(pmat.cell2cluster))>1){
    warning("Not all cells are randomly sampled. Consider increase the number of iterations specified in the 'subsamp.iter' argument.")
  }
  pmat.cell2cluster[is.na(pmat.cell2cluster)] <- 0

  ## all outputs
  rmax.cell2cluster <- apply(pmat.cell2cluster,1,max)
  match.cell2cluster <- clusterNames.ref[max.col(pmat.cell2cluster)]
  tab.query <- sapply(rstlst.query.pmat,nrow)
  query.cluster <- rep(names(tab.query),tab.query)
  output <- data.frame(query.cluster, match.cell2cluster, rmax.cell2cluster, pmat.cell2cluster, stringsAsFactors=FALSE) %>%
    rownames_to_column()
  colnames(output)[1:2] <- paste0("query.",c("cell","cluster"))

  return(output)
}


