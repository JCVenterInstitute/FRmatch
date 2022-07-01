
########################################################################################################
## filter_cluster()
########################################################################################################

filter_cluster <- function(sce.object, filter.size, filter.fscore=NULL, filter.nomarker=FALSE){

  ## cluster size
  tab <- table(colData(sce.object)$cluster_membership)

  ## filter by cluster size
  cluster.keep <- names(which(tab>=filter.size))

  ## filter by fscore
  if(!is.null(filter.fscore)){
    cluster.keep <- sce.object@metadata$f_score %>%
      filter(clusterName %in% cluster.keep & score>=filter.fscore) %>%
      pull(clusterName)
  }

  # ## colData
  # ind.keep <- colData(sce.object)$cluster_membership %in% cluster.keep
  # sce.object.filt <- sce.object[,ind.keep]
  # ## metaData
  # if(!is.null(sce.object@metadata$cluster_marker_info)){
  #   sce.object.filt@metadata$cluster_marker_info <- sce.object@metadata$cluster_marker_info %>%
  #     filter(clusterName %in% cluster.keep)
  # }
  # if(!is.null(sce.object@metadata$f_score)){
  #   sce.object.filt@metadata$f_score <- sce.object@metadata$f_score %>%
  #     filter(clusterName %in% cluster.keep)
  # }
  # if(!is.null(sce.object@metadata$cluster_order)){
  #   sce.object.filt@metadata$cluster_order <- base::intersect(sce.object@metadata$cluster_order, cluster.keep)
  # }
  # ## rowData
  # if(!is.null(sce.object@metadata$cluster_marker_info)){
  #   rowData(sce.object.filt)$NSF_markers <- rownames(sce.object.filt) %in% sce.object.filt@metadata$cluster_marker_info$markerGene
  # }

  sce.object.filt <- subset_by_cluster(sce.object, cluster.keep)

  ## output
  return(sce.object.filt)
}


########################################################################################################
## subset_by_cluster()
########################################################################################################

subset_by_cluster <- function(sce.object, cluster.keep){

  ## colData
  ind.keep <- colData(sce.object)$cluster_membership %in% cluster.keep
  sce.object.filt <- sce.object[,ind.keep]
  ## metaData
  if(!is.null(sce.object@metadata$cluster_marker_info)){
    sce.object.filt@metadata$cluster_marker_info <- sce.object@metadata$cluster_marker_info %>%
      filter(clusterName %in% cluster.keep)
  }
  if(!is.null(sce.object@metadata$f_score)){
    sce.object.filt@metadata$f_score <- sce.object@metadata$f_score %>%
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


########################################################################################################
## impute_zero()
########################################################################################################

# #' Imputation for zero expression of marker genes within cluster
# #'
# #' Impute the reference dataset only. Do NOT impute the query dataset. See 'Details' for imputation scheme.
# #'
# #' @param sce.object A \link[SingleCellExperiment]{SingleCellExperiment} data object for \code{FRmatch}.
# #' See example in \code{\link[FRmatch]{sce.example}}.
# #'
# #' @return A data object with imputed values for the zero expression of NS-Forest markers within their specific clusters.
# #'
# #' @details
# #' Imputation scheme:
# #' \itemize{
# #' \item Assume that the cluster-specific NS-Forest marker genes have non-zero expression in the cell type cluster that they characterize.
# #' \item For each marker gene, impute the zero values \emph{within} in the specific cluster.
# #' \item Draw single imputation value from the empirical distribution formed by non-zero values in the cluster.
# #' }
# #'
# #' @examples
# #' \dontrun{
# #' data("sce.example")
# #' plot_nonzero(sce.example)
# #' sce.example.imputed <- impute_zero(sce.example)
# #' plot_nonzero(sce.example.imputed)
# #' }


impute_zero <- function(sce.object){
  dat <- assay(sce.object)
  cluster_membership <- colData(sce.object)$cluster_membership
  cluster_marker_info <- metadata(sce.object)$cluster_marker_info
  cluster_names <- unique(cluster_marker_info$clusterName)
  cluster_marker_mat <- dat[cluster_marker_info$markerGene,]

  my.split <- function(cluster, markerGene, datmat, cluster_membership){
    cluster_names <- unique(cluster)
    datlst <- vector("list", length=length(cluster_names))
    names(datlst) <- cluster_names
    for (cl in cluster_names){
      datlst[[cl]] <- t(t(datmat)[cluster_membership==cl,markerGene[cluster==cl]]) #t(t()) to make sure if a vector is in row matrix
    }
    return(datlst)
  }
  out <- my.split(cluster_marker_info$clusterName, cluster_marker_info$markerGene, cluster_marker_mat, cluster_membership)

  ## imputation
  my.impute <- function(zz){
    ind.zero <- zz==0
    if(sum(ind.zero)>0){
      zz[ind.zero] <- pmax(0,rnorm(sum(ind.zero), mean(zz[!ind.zero]), ifelse(is.na(sd(zz[!ind.zero])), 0, sd(zz[!ind.zero]))))
    }
    else zz
    return(zz)
  }
  out.imputed <- lapply(out, function(z) t(apply(z,1,my.impute)))
  ## check
  # lapply(out.imputed, function(z) rowSums(z==0)) #check no more zeros

  ## merge back to the whole data matrix
  my.merge <- function(cluster, markerGene, datmat, cluster_membership, out.imputed){
    cluster_names <- names(out.imputed)
    for (cl in cluster_names){
      datmat[markerGene[cluster==cl],cluster_membership==cl] <- out.imputed[[cl]]
    }
    return(datmat)
  }
  dat.imputed <- my.merge(cluster_marker_info$clusterName, cluster_marker_info$markerGene, dat, cluster_membership, out.imputed)

  assay(sce.object) <- dat.imputed
  return(sce.object)
}
