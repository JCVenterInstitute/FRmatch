
#' FR-Match for single cell RNA-seq data
#'
#' This is a user-end wrapper function that implements the steps of cell type matching between two single cell RNA-seq experiments
#' (namely, \code{query} and \code{reference}) using the cluster-to-cluster FR-Match approach.
#'
#' @param sce.query Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for query experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param sce.ref Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for reference experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param use.cosine Boolean variable indicating if to use cosine distance in \code(\link[FRmatch]{FR-test}}. Default: \code{TRUE}.
# #' @param imputation INACTIVE. Boolean variable indicating if to impute expression zero values for the reference experiment. Default: \code{FALSE}.
# #' See details in \code{\link[FRmatch]{impute_dropout}}.
#' @param filter.size,filter.fscore,filter.nomarker Filtering out small/poor-quality/no-marker clusters. Default: \code{filter.size=5}, filter based on the number
#' of cells per cluster; \code{filter.fscore=NULL}, do not filter based on the F-beta score, otherwise specify a numeric value between 0 and 1;
#' \code{filter.nomarker=FALSE}, filter based on the boolean variable indicating if to filter reference clusters with no marker genes available in query.
#' @param add.pseudo.marker,pseudo.expr Adding pseudo marker to stabilize no expression clusters in the marker gene feature space.
#' Default: \code{add.pseudo.marker=FALSE}, boolean. Pseudo marker expression values are drawn from uniform distribution from 0 to \code{pseudo.expr}.
#' Default: \code{pseudo.expr=1}, numeric, for the min-max scaled data after normalization.
#' @param subsamp.size,subsamp.iter,subsamp.seed Numeric variables for iterative subsampling size, number of iterations, and random seed for iterations. YMMV.
#' @param numCores Number of cores for parallel computing. Default: \code{NULL}, use the maximum number of cores detected by \code{\link[parallel]{detectCores}}.
#' Otherwise, specify by an integer value.
#' @param prefix Prefix names for query and reference clusters. Default: \code{prefix=c("query.", "ref.")}.
#' @param verbose Numeric value indicating levels of details to be printed. Default: \code{1}, only print major steps.
#' If \code{0}, no verbose; if \code{2}, print all, including warnings.
#' @param return.all Boolean variable indicating if to return all results (such as runs, etc.). Default: \code{FALSE}.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FRtest}}, including \code{use.cosine}.
#'
#' @return A list of:
#' \item{settings}{Record of customized parameter settings specified in the function.}
#' \item{pmat}{A matrix of p-values. Rows are reference clusters, and columns are query clusters.}
#' \item{statmat}{A matrix of FR statistics. Rows are reference clusters, and columns are query clusters.}
#' If \code{return.all = TRUE}, more intermediate results are returned.
#'
#' @seealso Visualization of matching results using \code{\link[FRmatch]{plot_FRmatch}}, \code{\link[FRmatch]{plot_bi_FRmatch}}.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' FRmatch(sce.example, sce.example)
#' }
#'
#' @importFrom pbmcapply pbmcmapply
#'
#' @export

FRmatch <- function(sce.query, sce.ref, use.cosine=TRUE,  #imputation=FALSE,
                    filter.size=5, filter.fscore=NULL, filter.nomarker=FALSE, #filtering clusters
                    add.pseudo.marker=FALSE, pseudo.expr=1, #adding pseudo marker
                    subsamp.size=20, subsamp.iter=1000, subsamp.seed=1, #subsampling
                    numCores=NULL, prefix=c("query.", "ref."),
                    verbose=1, return.all=FALSE, ...){

  ## check data object
  if(verbose>0) cat("* Check query data object. \n")
  sce.query <- check_data_object(sce.query, is.ref=FALSE, verbose=(verbose>1))
  if(verbose>0) cat("* Check reference data object. \n")
  sce.ref <- check_data_object(sce.ref, verbose=(verbose>1))

  ## find common marker genes
  markergenes <- rownames(sce.ref)[rowData(sce.ref)$marker_gene==1]
  markergenes.common <- intersect(markergenes, rownames(sce.query))
  if(verbose>0) cat("* Feature selection:", length(markergenes.common), "out of", length(markergenes),
                    "reference marker genes are presented in the query experiment. \n")
  if(length(markergenes.common)!=length(markergenes)){
    if(verbose>1) warning(paste(setdiff(markergenes,markergenes.common),collapse = ", ")," not found.")
    # cat(paste(setdiff(markergenes,markergenes.common),collapse = ", "),"not found. \n")
  }

  ## filtering ref clusters without marker genes available in query
  if(filter.nomarker){
    if(verbose>0) cat("* Filtering reference clusters with no marker genes available in query. \n")
    ref.cluster.keep <- sce.ref@metadata$cluster_marker_info %>% filter(markerGene %in% markergenes.common) %>% pull(clusterName) %>% unique()
    ref.cluster.nomarker <- setdiff(unique(colData(sce.ref)$cluster_membership), ref.cluster.keep)
    sce.ref <- FRmatch:::subset_by_cluster(sce.ref, ref.cluster.keep)
    if(length(ref.cluster.nomarker)>0){
      if(verbose>1) warning(paste(ref.cluster.nomarker, collapse = ", ")," with no marker gene available in query.")
    }
  }

  ## filtering small or low fscore clusters
  if(verbose>0) cat("* Filtering small clusters: query and reference clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)){
    if(verbose>0) cat("* Filtering low F-beta score clusters: reference cluster with F-beta score <", filter.fscore, "are not considered. \n")
  }
  sce.query <- filter_cluster(sce.query, filter.size=filter.size) #only filter on size, not fscore
  sce.ref <- filter_cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)

  ## imputation
  # if(imputation){
  #   sce.ref <- impute_zero(sce.ref)
  #   cat("Imputation is applied. \n")
  # }

  ## get expr data and dimension reduction by selecting common marker genes
  querydat.reduced <- logcounts(sce.query)[markergenes.common,]
  refdat.reduced <- logcounts(sce.ref)[markergenes.common,]

  ## if to add pseudo marker to give signals in query clusters with almost no expression
  if(add.pseudo.marker){
    if(verbose>0) cat("* Adding pseudo marker to stablize no expression clusters.\n")
    querydat.reduced <- rbind(querydat.reduced, runif(ncol(querydat.reduced), 0, pseudo.expr))
    refdat.reduced <- rbind(refdat.reduced, runif(ncol(refdat.reduced), 0, pseudo.expr))
  }

  ## extract cluster info from sce.objects
  membership.query <- colData(sce.query)$cluster_membership
  membership.ref <- colData(sce.ref)$cluster_membership
  order.query <- sce.query@metadata$cluster_order
  order.ref <- sce.ref@metadata$cluster_order

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

  if(verbose>0) cat("** method = cluster2cluster", "| subsamp.size =", subsamp.size, "| subsamp.iter =", subsamp.iter, "\n")

  results <- pbmcmapply(
    function(samp1,samp2){
      set.seed(subsamp.seed)
      FRtest_subsamp(samp1, samp2, use.cosine=use.cosine, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, ...)
    },
    paired.datlst.query, paired.datlst.ref,
    mc.cores = numCores)

  if(verbose>0) cat("* Done parallel and returning results. \n")

  ## outputs
  pmat <- matrix(results["p.value",], nrow=ncluster.ref)
  statmat <- matrix(results["stat",], nrow=ncluster.ref)
  rownames(statmat) <- rownames(pmat) <- clusterNames.ref
  colnames(statmat) <- colnames(pmat) <- clusterNames.query

  ## return
  settings <- list(filter.size=filter.size, filter.fscore=filter.fscore, filter.nomarker=filter.nomarker,
                   add.pseudo.marker=add.pseudo.marker, pseudo.expr=pseudo.expr,
                   use.cosine=use.cosine, method="cluster2cluster",
                   subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, subsamp.seed=subsamp.seed)
  if(return.all){
    ## all results
    all.results <- as.data.frame(t(results)) %>%
      mutate(query.cluster=rep(clusterNames.query, each=ncluster.ref),
             ref.cluster=rep(clusterNames.ref, times=ncluster.query)) %>%
      select(query.cluster, ref.cluster, p.value, stat, runs, runs.query=runs.samp1, runs.ref=runs.samp2)
    return(list("settings"=settings, "pmat"=pmat, "statmat"=statmat, "all.results"=all.results))
  } else {
    return(list("settings"=settings, "pmat"=pmat, "statmat"=statmat))
  }
}


