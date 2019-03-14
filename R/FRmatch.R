#######################
## TO DO:
## 1. bootstrap
## 2. multiple hypothesis tesing correction
#######################

#' Cross-comparison of two single-cell RNA-seq experiments using Friedman-Rafsky test
#'
#' This is a user-end function that wraps the steps to cross-compare two single-cell RNA-seq experiments (namely, query and reference)
#' with cluster information and NS-Forest marker genes using Friedman-Rafsky (FR) test based methods.
#' This function inputs two objects of \link[SingleCellExperiment]{SingleCellExperiment} class with necessary information,
#' and returns FR statistics (a.k.a. scores) and p-values of the FR test.
#'
#' @param sce.query Query experiment object of \code{SingleCellExperiment} class with necessary information.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param sce.ref Reference experiment object of \code{SingleCellExperiment} class with necessary information.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param method RIGHT NOW ONLY IMPLEMENTED THE STANDARD FR TEST.
#' @param imputation Logical variable indicating if to impute dropout values in the reference experiment. Default: \code{FALSE}.
#' See details in \code{\link[FRmatch]{impute_dropout}}.
#' @param numCores \code{NULL} or an integer that specifies the number of cores to be used in parallel computing.
#' Default: \code{NULL}, which sets \code{numCores = detectCores()}. See more in \code{\link[parallel]{detectCores}}.
#' @param return.all Logical variable indicating if to return all results (such as number of subtrees, etc.).
#' Default: \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FR.test}}.
#'
#' @return A list of FR test results:
#' \item{statmat}{A matrix of all FR statistics for pairwise comparison. Rows are reference clusters, and columns are query clusters.}
#' \item{pmat}{A matrix of all p-values for pairwise comparison. Rows are reference clusters, and columns are query clusters.}
#' If \code{return.all = TRUE}, more results are returned in \code{all.restuls} element of the list.
#'
#' @author Yun Zhang, \email{zhangy@jcvi.org}; Brian Aevermann, \email{baeverma@jcvi.org}; Richard Scheuermann, \email{RScheuermann@jcvi.org}.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.
#'
#' @examples
#' @export

FRmatch <- function(sce.query, sce.ref, method="none", imputation=FALSE, numCores=NULL, return.all=FALSE, ...){
  ## check data object
  cat("Check query data object: ")
  sce.query <- check_query_data(sce.query)
  cat("Check reference data object: ")
  sce.ref <- check_data(sce.ref)

  ## imputation
  if(imputation){
    sce.ref <- impute_dropout(sce.ref)
    cat("Imputation done. \n")
  }

  ## extract info from sce.objects
  querydat <- logcounts(sce.query) #matrix
  refdat <- logcounts(sce.ref)
  membership.query <- colData(sce.query)$cluster_membership
  membership.ref <- colData(sce.ref)$cluster_membership
  order.query <- metadata(sce.query)$cluster_order
  order.ref <- metadata(sce.ref)$cluster_order

  ## dimension reduction by selecting only marker genes
  markergenes <- metadata(sce.ref)$cluster_marker_info$markerGene
  markergenes.common <- intersect(markergenes, rownames(sce.query))
  querydat.reduced <- querydat[markergenes.common,]
  refdat.reduced <- refdat[markergenes.common,]
  cat("Dimension reduction done:", length(markergenes.common), "out of", length(markergenes),
      "marker genes in the reference experiment are presented in the query experiment. \n")
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

  # if(method=="none"){
  #   results <- lapply(1:length(datlst.query),
  #                     function(i) sapply(1:length(datlst.ref), function(j){
  #                       # FR.test(datlst.query[[i]], datlst.ref[[j]], ...)
  #                       FR.test(datlst.query[[i]], datlst.ref[[j]])}
  #                     )
  #   )
  # }
  # # if(method=="boot"){
  # #   results <- lapply(1:length(lst.query),
  # #                     function(i) sapply(1:length(lst.ref),
  # #                                        function(j) FR.test.boot(lst.query[[i]], lst.ref[[j]], ...)
  # #                     )
  # #   )
  # # }
  # ## assign names
  # names(results) <- names(datlst.query)
  # restuls <- lapply(results, function(zz){
  #   colnames(zz) <- names(datlst.ref)
  #   zz})
  # ## outputs
  # statmat <- sapply(restuls, "[", "stat",)
  # pmat <- sapply(restuls, "[", "p.value",)
  # # colnames(statmat) <- colnames(pmat) <- colnames(plst)
  # # rownames(statmat) <- rownames(pmat) <- rownames(plst)
  # return(list("statmat"=statmat, "pmat"=pmat))

  ##------ USE PARALLEL COMPUTING!!! ------##

  ## use all available cores if not specified
  if(is.null(numCores)) numCores <- detectCores()

  ## prepare data for each combination of pairs between query and ref clusters
  paired.datlst.query <- rep(datlst.query, each=ncluster.ref)
  paired.datlst.ref <- rep(datlst.ref, times=ncluster.query)

  if(method=="none"){
    results <- mcmapply(function(samp1,samp2){FR.test(samp1, samp2, ...)},
                        paired.datlst.query, paired.datlst.ref,
                        mc.cores = numCores)
  }

  ##------------##
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
    return(list("statmat"=statmat, "pmat"=pmat, "all.results"=all.results))
  }
  else return(list("statmat"=statmat, "pmat"=pmat))
}


