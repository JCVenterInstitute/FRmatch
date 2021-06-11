
#' Normalization function for FR-Match pipeline
#'
#' @param sce.object Data object.
#' @param scale Boolean variable indicating if to scale the per gene expression range to 0 and 1. Default: \code{TRUE}.
#' @param norm.by If \code{norm.by=NULL} (default), do nothing. If \code{norm.by="median"}, to weight each gene by its median then rescale to unit range.
#' If \code{norm.by="mean"}, to weight each gene by its mean then rescale to unit range.
#'
#' @return A data object with normalized expression data.
#'
#' @export

normalization <- function(sce.object, scale=TRUE, norm.by=NULL){

  ## scale to 0 and 1
  if(scale){
    mat <- assay(sce.object)
    mat.scale <- sweep(mat,1,rowMax(mat),"/")
    assay(sce.object) <- mat.scale
  }

  if(is.null(norm.by)) norm.by="none"
  ## normalize by row median
  if(norm.by=="median"){
    mat.norm <- assay(sce.object)
    clusters <- unique(colData(sce.object)$cluster_membership)
    for(cl in clusters){
      ind.col <- colData(sce.object)$cluster_membership==cl
      mat.cl <- assay(sce.object)[,ind.col] %>%
        as.matrix() #in case there is only one cell in the cluster
      mat.cl.norm <- sweep(mat.cl,1,rowMedians(mat.cl),"*")
      ## final rescaling to 0 and 1
      mat.cl.norm.rescale <- sweep(mat.cl.norm,1,max(mat.cl.norm),"/")
      mat.norm[,ind.col] <- mat.cl.norm.rescale
      assay(sce.object) <- mat.norm
    }
  }
  ## normalize by row mean
  if(norm.by=="mean"){
    mat.norm <- assay(sce.object)
    clusters <- unique(colData(sce.object)$cluster_membership)
    for(cl in clusters){
      ind.col <- colData(sce.object)$cluster_membership==cl
      mat.cl <- assay(sce.object)[,ind.col] %>%
        as.matrix() #in case there is only one cell in the cluster
      mat.cl.norm <- sweep(mat.cl,1,rowMeans(mat.cl),"*")
      ## final rescaling to 0 and 1
      mat.cl.norm.rescale <- sweep(mat.cl.norm,1,max(mat.cl.norm),"/")
      mat.norm[,ind.col] <- mat.cl.norm.rescale
      assay(sce.object) <- mat.norm
    }
  }

  return(sce.object)
}

