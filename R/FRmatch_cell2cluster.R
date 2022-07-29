
#' FR-Match cell-to-cluster matching
#'
#' This is a user-end wrapper function that implements the steps of cell type matching between two single cell RNA-seq experiments
#' (namely, \code{query} and \code{reference}) using the cell-to-cluster FR-Match approach that matches each query cell to a reference cluster.
#'
#' @param sce.query Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for query experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param sce.ref Data object of the \link[SingleCellExperiment]{SingleCellExperiment} class for reference experiment.
#' See details in \code{\link[FRmatch]{sce.example}}.
# #' @param imputation INACTIVE. Boolean variable indicating if to impute expression zero values for the reference experiment. Default: \code{FALSE}.
# #' See details in \code{\link[FRmatch]{impute_dropout}}.
#' @param feature.selection Which set of features to use for the matching space? Default: \code{feature.selection="reference.markers"}, use reference marker genes.
#' If \code{feature.selection="query.genes"}, use all query genes as the feature space, e.g. the query genes are probe genes for spatial transcriptomics experiment
#' @param filter.size,filter.fscore,filter.nomarker Filtering out small/poor-quality/no-marker clusters. Default: \code{filter.size=5}, filter based on the number
#' of cells per cluster; \code{filter.fscore=NULL}, do not filter based on the F-beta score, otherwise specify a numeric value between 0 and 1;
#' \code{filter.nomarker=FALSE}, filter based on the boolean variable indicating if to filter reference clusters with no marker genes available in query
#' in the case \code{feature.selection="reference.markers"}.
#' @param add.pseudo.marker,pseudo.expr Adding pseudo marker to stabilize no expression clusters in the marker gene feature space.
#' Default: \code{add.pseudo.marker=FALSE}, boolean. Pseudo marker expression values are drawn from uniform distribution from 0 to \code{pseudo.expr}.
#' Default: \code{pseudo.expr=1}, numeric, for the min-max scaled data after normalization.
#' @param subsamp.size,subsamp.iter,subsamp.seed Numeric variables for iterative subsampling size, number of iterations, and random seed for iterations. YMMV.
#' @param numCores Number of cores for parallel computing. Default: \code{NULL}, use the maximum number of cores detected by \code{\link[parallel]{detectCores}}.
#' Otherwise, specify by an integer value.
#' @param prefix Prefix names for query and reference clusters. Default: \code{prefix=c("query.", "ref.")}.
#' @param verbose Numeric value indicating levels of details to be printed. Default: \code{1}, only print major steps.
#' If \code{0}, no verbose; if \code{2}, print all, including warnings.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FRtest}}, including \code{use.cosine}.
#'
#' @return A list of:
#' \item{settings}{Record of customized parameter settings specified in the function.}
#' \item{pmat}{A cell-by-cluster (a.k.a. query cell by reference cluster) matrix of p-values retained from the iterative procedure.}
#' \item{cell2cluster}{A data frame of cell-to-cluster matches summarized from the \code{pmat}.}
#' Columns in \code{cell2cluster} are:
#' \item{query.cell}{Query cell ID.}
#' \item{query.cluster}{Cluster membership of query cells.}
#' \item{match}{Matched reference cluster for the query cell.}
#' \item{score}{Confidence score of \code{match.cell2cluster}, which is the maximum value of the corresponding row in \code{pmat}.}
#'
#' @details
#' This implementation is \code{FRmatch} with an iterative subsampling scheme, which is a bootstrap-like approach to randomly select a smaller
#' set of cells from a query cluster and quantify the confidence score of the selected cells belonging to certain reference cell type
#' using the p-value outputted from \code{FRtest} (i.e. a larger p-value indicates higher probability of a match, and vice versa).
#'
#' This function assigns the cluster-level p-value to each selected query cell, and updates the assigned p-value if the query cell is reselected
#' from the iterative procedure and assigned a higher p-value. The output from this implementation includes a cell-by-cluster
#' (a.k.a. query cell by reference cluster) matrix of p-values.
#'
#' @seealso Visualization of matching results using \code{\link[FRmatch]{plot_FRmatch_cell2cluster}}}.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' FRmatch_cell2cluster(sce.example, sce.example)
#' }
#'
#' @export

FRmatch_cell2cluster <- function(sce.query, sce.ref, use.cosine=TRUE,  #imputation=FALSE,
                                 feature.selection="reference.markers", #feature selection
                                 filter.size=5, filter.fscore=NULL, filter.nomarker=FALSE, #filtering clusters
                                 add.pseudo.marker=FALSE, pseudo.expr=1, #adding pseudo marker
                                 subsamp.size=10, subsamp.iter=2000, subsamp.seed=1, #subsampling
                                 numCores=NULL, prefix=c("query.", "ref."),
                                 verbose=1, ...){

  #######################
  ## check data object ##
  #######################
  if(verbose>0) cat("* Check query data object. \n")
  sce.query <- check_data_object(sce.query, verbose=(verbose>1))
  if(verbose>0) cat("* Check reference data object. \n")
  sce.ref <- check_data_object(sce.ref, verbose=(verbose>1))

  #######################
  ## feature selection ##
  #######################
  ## use query gene space
  if(feature.selection=="query.genes"){
    markergenes <- rownames(sce.query)
    markergenes.common <- intersect(markergenes, rownames(sce.ref))
    if(verbose>0) cat("* Feature selection: gene space of the query experiment. \n",
                      "**", length(markergenes.common), "out of", length(markergenes), "query genes are common in the reference experiment. \n")
    if(length(markergenes.common)!=length(markergenes)){
      if(verbose>1) warning(paste(setdiff(markergenes,markergenes.common),collapse = ", ")," not found.")
      # cat(paste(setdiff(markergenes,markergenes.common),collapse = ", "),"not found. \n")
    }
  }

  ## use reference markers
  if(feature.selection=="reference.markers"){
    ## find common marker genes
    markergenes <- rownames(sce.ref)[rowData(sce.ref)$marker_gene==1]
    markergenes.common <- intersect(markergenes, rownames(sce.query))
    if(verbose>0) cat("* Feature selection: marker genes of the reference experiment. \n",
                      "**", length(markergenes.common), "out of", length(markergenes), "reference marker genes are presented in the query experiment. \n")
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
  }

  #######################
  ## cluster filtering ##
  #######################
  ## filtering small or low fscore clusters
  if(verbose>0) cat("* Filtering small clusters: query and reference clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)){
    if(verbose>0) cat("* Filtering low F-beta score clusters: reference cluster with F-beta score <", filter.fscore, "are not considered. \n")
  }
  sce.query <- filter_cluster(sce.query, filter.size=filter.size) #only filter on size, not fscore
  sce.ref <- filter_cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)

  ################
  ## imputation ##
  ################
  # if(imputation){
  #   sce.ref <- impute_zero(sce.ref)
  #   cat("Imputation is applied. \n")
  # }

  ##################
  ## reduced data ##
  ##################
  ## get expr data and dimension reduction by selecting common marker genes
  querydat.reduced <- logcounts(sce.query)[markergenes.common,]
  refdat.reduced <- logcounts(sce.ref)[markergenes.common,]

  #######################
  ## add pseudo marker ##
  #######################
  ## if to add pseudo marker to give signals in query clusters with almost no expression
  if(add.pseudo.marker){
    if(verbose>0) cat("* Adding pseudo marker to stablize no expression clusters.\n")
    querydat.reduced <- rbind(querydat.reduced, runif(ncol(querydat.reduced), 0, pseudo.expr))
    refdat.reduced <- rbind(refdat.reduced, runif(ncol(refdat.reduced), 0, pseudo.expr))
  }

  ###############
  ## data prep ##
  ###############
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
  if(verbose>0) cat(" ** Parallel computing with", numCores, "cores. \n")

  ## prepare data for each combination of query cluster and ref cluster pair
  paired.datlst.query <- rep(datlst.query, each=ncluster.ref)
  paired.datlst.ref <- rep(datlst.ref, times=ncluster.query)
  ##---------------------------------------##

  if(verbose>0) cat(" ** method = cell2cluster", "| subsamp.size =", subsamp.size, "| subsamp.iter =", subsamp.iter, "\n")

  results <- pbmcmapply(
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
  score.cell2cluster <- apply(pmat.cell2cluster,1,max)
  match.cell2cluster <- clusterNames.ref[max.col(pmat.cell2cluster)]
  query.cell <- rownames(pmat.cell2cluster)
  ncell.query <- sapply(rstlst.query.pmat,nrow)
  query.cluster <- rep(names(ncell.query),ncell.query)
  output <- data.frame(query.cell, query.cluster, "match"=match.cell2cluster, "score"=score.cell2cluster, stringsAsFactors=FALSE)

  ## return
  settings <- list(filter.size=filter.size, filter.fscore=filter.fscore, filter.nomarker=filter.nomarker,
                   add.pseudo.marker=add.pseudo.marker, pseudo.expr=pseudo.expr,
                   use.cosine=use.cosine, method="cluster2cluster",
                   subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, subsamp.seed=subsamp.seed)
  return(list("settings"=settings, "pmat"=pmat.cell2cluster, "cell2cluster"=output))
}


