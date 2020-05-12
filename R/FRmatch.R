
#' Cluster-to-cluster cell type matching algorithm for single-cell RNA-seq data
#'
#' This is a user-end function that wraps up the steps of matching cell type clusters between two single-cell RNA-seq experiments
#' (namely, \code{query} and \code{reference}) with \emph{clustered} expression data and \emph{informative} marker genes
#' using the Friedman-Rafsky (FR) statistical test.
#' This function inputs two data objects of the \link[SingleCellExperiment]{SingleCellExperiment} class,
#' and outputs the FR statistics, p-values, and optionally, all intermediate restuls from the FR test.
#'
#' @param sce.query Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for query experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param sce.ref Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for reference experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param imputation INACTIVE. Logical variable indicating if to impute expression zero values for the reference experiment. Default: \code{FALSE}.
#' See details in \code{\link[FRmatch]{impute_dropout}}.
#' @param filter.size,filter.fscore Filtering small/poor-quality clusters. Default: \code{filter.size=10}, filter based on the number
#' of cells per cluster; \code{filter.fscore=NULL}, filter based on the f-score associated with the cell cluster if available (numeric value).
#' @param method Methods for the FR test. Default: \code{method="subsampling"} is to iteratively subsample equal number of cells (i.e. cluster size)
#' from the query and reference clusters, and then perform the FR test. Option: \code{method="none"} is the FR test with no modification.
#' @param subsamp.size,subsamp.iter,subsamp.seed Cluster size, number of iterations, and random seed for \code{method="subsampling"}.
#' Default: \code{10, 1001, 1}, respectively.
#' @param numCores Number of cores for parallel computing.
#' Default: \code{NULL}, use the maximum number of cores detected by \code{\link[parallel]{detectCores}} if not specified (an integer).
#' @param prefix Prefix names for query and reference clusters. Default: \code{prefix=c("query.", "ref.")}.
#' @param verbose Numeric value indicating levels of details to be printed. Default: \code{1}, only print major steps.
#' If \code{0}, no verbose; if \code{2}, print all.
#' @param return.all Logical variable indicating if to return all results (such as runs, etc.). Default: \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FR.test}}.
#'
#' @return A list of:
#' \item{parameters}{Call of key parameters used in the analysis.}
#' \item{pmat}{A matrix of p-values. Rows are reference clusters, and columns are query clusters.}
#' \item{statmat}{A matrix of FR statistics. Rows are reference clusters, and columns are query clusters.}
#' If \code{return.all = TRUE}, more intermediate results are returned.
#'
#' @author Yun Zhang, \email{zhangy@jcvi.org};
#' Brian Aevermann, \email{baeverma@jcvi.org};
#' Richard Scheuermann, \email{RScheuermann@jcvi.org}.
#'
#' @seealso Visualization of matching results using \code{\link[FRmatch]{plot_FRmatch}}, \code{\link[FRmatch]{plot_bilateral_FRmatch}}.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' FRmatch(sce.example, sce.example)
#' }
#'
#' @export

FRmatch <- function(sce.query, sce.ref, #imputation=FALSE,
                    filter.size=10, filter.fscore=NULL, #filtering clusters
                    method="subsampling", subsamp.size=10, subsamp.iter=1001, subsamp.seed=1, #subsampling
                    numCores=NULL, prefix=c("query.", "ref."),
                    verbose=1, return.all=FALSE, ...){

  ## check data object
  if(verbose>0) cat("* Check query data object. \n")
  sce.query <- check_data_object(sce.query, verbose=(verbose>1))
  if(verbose>0) cat("* Check reference data object. \n")
  sce.ref <- check_data_object(sce.ref, verbose=(verbose>1))

  ## filtering small or low fscore clusters
  if(verbose>0) cat("* Filtering small clusters: query and reference clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)){
    if(verbose>0) cat("* Filtering low F-score clusters: reference cluster with f-score <", filter.fscore, " are not considered. \n")
  }
  sce.query <- filter.cluster(sce.query, filter.size=filter.size, filter.fscore=NULL)
  sce.ref <- filter.cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)

  ## imputation
  # if(imputation){
  #   sce.ref <- impute.zero(sce.ref)
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
  ##------------##

  if(method=="none"){
    if(verbose>0) cat("** method =", method, "\n")
    results <- mcmapply(function(samp1,samp2){FR.test(samp1, samp2, ...)},
                        paired.datlst.query, paired.datlst.ref,
                        mc.cores = numCores)
  }
  if(method=="subsampling"){
    if(verbose>0) cat("** method =", method, "| subsamp.size =", subsamp.size, "| subsamp.iter =", subsamp.iter, "\n")
    set.seed(subsamp.seed)
    results <- mcmapply(
      function(samp1,samp2){
        FR.test.subsamp(samp1, samp2, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, ...)
      },
      paired.datlst.query, paired.datlst.ref,
      mc.cores = numCores)
  }
  if(verbose>0) cat("* Done parallel and returning results. \n")

  ## outputs
  pmat <- matrix(results["p.value",], nrow=ncluster.ref)
  statmat <- matrix(results["stat",], nrow=ncluster.ref)
  rownames(statmat) <- rownames(pmat) <- clusterNames.ref
  colnames(statmat) <- colnames(pmat) <- clusterNames.query

  ## return
  parameters <- list(filter.size=filter.size, filter.fscore=filter.fscore, method=method, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter)
  if(return.all){
    ## all results
    all.results <- as.data.frame(t(results)) %>%
      mutate(query.cluster=rep(clusterNames.query, each=ncluster.ref),
             ref.cluster=rep(clusterNames.ref, times=ncluster.query)) %>%
      select(query.cluster, ref.cluster, p.value, stat, runs, runs.query=runs.samp1, runs.ref=runs.samp2)
    return(list("parameters"=parameters, "pmat"=pmat, "statmat"=statmat, "all.results"=all.results))
  } else {
    return(list("parameters"=parameters, "pmat"=pmat, "statmat"=statmat))
  }
}


