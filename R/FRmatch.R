
#' Cell type cluster matching of single-cell RNA-seq data using Friedman-Rafsky test
#'
#' This is a user-end function that wraps the steps of matching cell type clusters between two single-cell RNA-seq experiments
#' (namely, \code{query} and \code{reference}) with clustered expression data and NS-Forest marker genes using Friedman-Rafsky (FR) test.
#' This function inputs two data objects of the \link[SingleCellExperiment]{SingleCellExperiment} data class,
#' and outputs the FR statistics, p-values, and optional intermediated restuls from the FR test.
#'
#' @param sce.query Data object of the \link[SingleCellExperiment]{SingleCellExperiment} data class for query experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param sce.ref A Data object of the \link[SingleCellExperiment]{SingleCellExperiment} data class for reference experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param imputation Logical variable indicating if to impute expression zero values for the reference experiment. Default: \code{FALSE}.
#' See details in \code{\link[FRmatch]{impute_dropout}}.
#' @param filter.size,filter.fscore Criteria for filtering small/poor clusters. Default: \code{filter.size=10}, filter based on the number
#' of cells per cluster; \code{filter.fscore=NULL}, filter based on the f-measure associated with NS-Forest markers (\code{NULL} or numeric value).
#' @param method Adapted methods of FR test. Default: \code{method="subsampling"} is to apply equal size subsampling of the comparing clusters for
#' the FR test. Option: \code{method="none"} is the standard FR test without modification.
#' @param subsamp.size,subsamp.iter,subsamp.seed Cluster size, number of iterations, and random seed for subsampling.
#' Default: \code{10, 1001, 1}, respectively.
#' @param numCores \code{NULL} or an integer that specifies the number of cores for parallel computing.
#' Default: \code{NULL}, is to set the maximum number of cores detected by \code{\link[parallel]{detectCores}} if not specified.
#' @param return.all Logical variable indicating if to return all results (such as number of runs, etc.). Default: \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FR.test}}.
#'
#' @return A list of FR test results:
#' \item{pmat}{A matrix of p-values for pairwise comparisons. Rows are reference clusters, and columns are query clusters.}
#' \item{statmat}{A matrix of FR statistics for pairwise comparisons. Rows are reference clusters, and columns are query clusters.}
#' If \code{return.all = TRUE}, more FR test intermediate results are returned.
#'
#' @author Yun Zhang, \email{zhangy@jcvi.org};
#' Brian Aevermann, \email{baeverma@jcvi.org};
#' Richard Scheuermann, \email{RScheuermann@jcvi.org}.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} data class.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' FRmatch(sce.example, sce.example)
#' }
#' @export

FRmatch <- function(sce.query, sce.ref, imputation=FALSE,
                    filter.size=10, filter.fscore=NULL, #filtering clusters
                    method="subsampling", subsamp.size=10, subsamp.iter=1001, subsamp.seed=1, #subsampling
                    numCores=NULL, return.all=FALSE, ...){

  ## check data object
  cat("Check query data object: ")
  sce.query <- check_query_data(sce.query)
  cat("Check reference data object: ")
  sce.ref <- check_data(sce.ref)

  ## filtering small or low fscore clusters
  cat("Filtering small clusters: clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)) cat("Filtering low f-measure clusters in the reference experiment only. \n")
  sce.ref <- filter.cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)
  sce.query <- filter.cluster(sce.query, filter.size=filter.size, filter.fscore=NULL)

  ## imputation
  if(imputation){
    sce.ref <- impute.zero(sce.ref)
    cat("Imputation done. \n")
  }

  ## extract info from sce.objects
  querydat <- logcounts(sce.query) #matrix
  refdat <- logcounts(sce.ref)
  membership.query <- SingleCellExperiment::colData(sce.query)$cluster_membership
  membership.ref <- SingleCellExperiment::colData(sce.ref)$cluster_membership
  order.query <- sce.query@metadata$cluster_order
  order.ref <- sce.ref@metadata$cluster_order

  ## dimension reduction by selecting only marker genes
  markergenes <- unique(metadata(sce.ref)$cluster_marker_info$markerGene)
  markergenes.common <- intersect(markergenes, rownames(sce.query))
  querydat.reduced <- querydat[markergenes.common,]
  refdat.reduced <- refdat[markergenes.common,]
  cat("Dimension reduction done:", length(markergenes.common), "out of", length(markergenes),
      "unique marker genes in the reference experiment are presented in the query experiment. \n")
  if(length(markergenes.common)!=length(markergenes)){
    warning(paste(setdiff(markergenes,markergenes.common),collapse = ", ")," not found in sce.query.")
  }

  ## split data by cluster membership and store in a list
  datlst.query <- lapply(split(as.data.frame(t(querydat.reduced)), membership.query),t)
  datlst.query <- datlst.query[order.query]
  clusterNames.query <- paste0("query.",names(datlst.query))
  names(datlst.query) <- clusterNames.query
  datlst.ref <- lapply(split(as.data.frame(t(refdat.reduced)), membership.ref),t)
  datlst.ref <- datlst.ref[order.ref]
  clusterNames.ref <- paste0("ref.",names(datlst.ref))
  names(datlst.ref) <- clusterNames.ref

  ncluster.ref <- length(datlst.ref)
  ncluster.query <- length(datlst.query)
  cat("Comparing", ncluster.query, "query clusters with", ncluster.ref,
      "reference clusters...")

  ###########################################
  ## FR comparison between two experiments ##
  ###########################################

  ##------ USE PARALLEL COMPUTING!!! ------##
  ## use all available cores if not specified
  if(is.null(numCores)) numCores <- detectCores()

  ## prepare data for each combination of pairs between query and ref clusters
  paired.datlst.query <- rep(datlst.query, each=ncluster.ref)
  paired.datlst.ref <- rep(datlst.ref, times=ncluster.query)
  ##------------##

  if(method=="none"){
    results <- mcmapply(function(samp1,samp2){FR.test(samp1, samp2, ...)},
                        paired.datlst.query, paired.datlst.ref,
                        mc.cores = numCores)
  }
  if(method=="subsampling"){
    set.seed(subsamp.seed)
    results <- mcmapply(
      function(samp1,samp2){
        FRmatch:::FR.test.subsamp(samp1, samp2, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter)
      },
      paired.datlst.query, paired.datlst.ref,
      mc.cores = numCores)
  }
  cat("Done. \n")

  ## outputs
  statmat <- matrix(results["stat",], nrow=ncluster.ref)
  pmat <- matrix(results["p.value",], nrow=ncluster.ref)
  rownames(statmat) <- rownames(pmat) <- clusterNames.ref
  colnames(statmat) <- colnames(pmat) <- clusterNames.query

  ## return
  if(return.all){
    ## all results
    all.results <- as.data.frame(t(results)) %>%
      dplyr::mutate(query.cluster=rep(clusterNames.query, each=ncluster.ref),
                    ref.cluster=rep(clusterNames.ref, times=ncluster.query)) %>%
      dplyr::select(query.cluster, ref.cluster, p.value, stat, runs, runs.query=runs.samp1, runs.ref=runs.samp2)
    return(list("pmat"=pmat, "statmat"=statmat, "all.results"=all.results))
  }
  else return(list("pmat"=pmat, "statmat"=statmat))
}


