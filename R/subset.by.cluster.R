subset.by.cluster <- function(sce.object, cluster.keep){
  
  ## colData
  ind.keep <- colData(sce.object)$cluster_membership %in% cluster.keep
  sce.object.filt <- sce.object[,ind.keep]
  ## metaData
  if(!is.null(sce.object@metadata$cluster_marker_info)){
    sce.object.filt@metadata$cluster_marker_info <- sce.object@metadata$cluster_marker_info %>%
      filter(cluster %in% cluster.keep)
  }
  if(!is.null(sce.object@metadata$fscores)){
    sce.object.filt@metadata$fscores <- sce.object@metadata$fscores %>%
      filter(clusterName %in% cluster.keep)
  }
  if(!is.null(sce.object@metadata$cluster_order)){
    sce.object.filt@metadata$cluster_order <- base::intersect(sce.object@metadata$cluster_order, cluster.keep)
  }
  ## rowData
  if(!is.null(sce.object@metadata$cluster_marker_info)){
    rowData(sce.object.filt)$NSF_markers <- rownames(sce.object.filt) %in% sce.object.filt@metadata$cluster_marker_info$markerGene
  }
  
  ## output
  return(sce.object.filt)
}