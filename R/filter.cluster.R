
filter.cluster <- function(sce.object, filter.size=10, filter.fscore=NULL){
  tab <- table(colData(sce.object)$cluster_membership)
  ## filter by cluster size
  cluster.keep <- names(which(tab>=filter.size))
  ## filter by fscore
  if(!is.null(filter.fscore)){
    cluster.keep <- sce.object@metadata$fscores %>%
      filter(clusterName %in% cluster.keep & `f-measure`>=filter.fscore) %>%
      pull(clusterName)
  }

  ## colData
  ind.keep <- colData(sce.object)$cluster_membership %in% cluster.keep
  sce.object.filt <- sce.object[,ind.keep]
  ## metaData
  sce.object.filt@metadata$cluster_marker_info <- sce.object@metadata$cluster_marker_info %>%
    filter(cluster %in% cluster.keep)
  if(!is.null(sce.object@metadata$fscores)){
    sce.object.filt@metadata$fscores <- sce.object@metadata$fscores %>%
      filter(clusterName %in% cluster.keep)
  }
  sce.object.filt@metadata$cluster_order <- base::intersect(sce.object@metadata$cluster_order, cluster.keep)
  ## rowData
  rowData(sce.object.filt)$NSF_markers <- rownames(sce.object.filt) %in% sce.object.filt@metadata$cluster_marker_info$markerGene

  return(sce.object.filt)
}
