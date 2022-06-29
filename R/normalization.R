
#' Normalization function for FR-Match analysis
#'
#' @param sce.object Data object.
#' @param scale Boolean variable indicating if to scale the per gene expression to the range [0,1]. Default: \code{TRUE}.
#' @param norm.by If \code{norm.by=NULL} (default), do nothing. If \code{norm.by="median"}, weight each gene by its median.
#' If \code{norm.by="mean"}, weight each gene by its mean. After normalization, all values are re-scale to [0,1] range.
#'
#' @return A data object with normalized expression data.
#'
#' @export

normalization <- function(sce.object, scale=TRUE, norm.by=NULL){

  ## scale to 0 and 1 PER GENE
  if(scale){
    mat <- assay(sce.object)
    mat.scale <- sweep(mat,1,rowMax(mat),"/")
    assay(sce.object) <- mat.scale
  }

  ## assign a character value for the if() statements below
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
      ## final re-scaling to 0 and 1
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
      ## final re-scaling to 0 and 1
      mat.cl.norm.rescale <- sweep(mat.cl.norm,1,max(mat.cl.norm),"/")
      mat.norm[,ind.col] <- mat.cl.norm.rescale
      assay(sce.object) <- mat.norm
    }
  }

  return(sce.object)
}

