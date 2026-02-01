
#' Normalization function for FR-Match analysis.
#'
#' @param sce.object Data object.
#' @param scale Boolean variable indicating if to scale the PER GENE expression to the range [0,1]. Defualt: \code{TRUE}.
#' @param norm.by Normalization options PER CLUSTER. Choose from \code{"softmax", "mean", "none"}.
#' If \code{"softmax"} (default), use the \code{\link[mclust]{softmax}} function to normalize all values globally. Afterwards, each cell is rescaled to unit length (i.e., L2 normalization).
#' If \code{"mean"}, each gene is weighted to de-noise lowly expressed genes. Afterwards, all values are rescaled to [0,1].
#' If \code{"none"}, do nothing.
#'
#' @return A data object with normalized expression data.
#'
#' @importFrom mclust softmax
#'
#' @export

normalization <- function(sce.object, scale=TRUE, norm.by="mean", impute.zero=FALSE){

  ## scale to 0 and 1 PER GENE globally
  if(scale){
    mat <- assay(sce.object)
    mat.scale <- sweep(mat,1,rowMax(mat),"/")
    mat.scale[is.na(mat.scale)] <- 0 #if rowMax(mat)==0
    assay(sce.object) <- mat.scale
  }

  ## normalize by softmax PER CLUSTER
  if(norm.by=="softmax"){
    mat.norm <- assay(sce.object)
    clusters <- unique(colData(sce.object)$cluster_membership)
    for(cl in clusters){
      ind.col <- colData(sce.object)$cluster_membership==cl
      mat.cl <- mat.norm[,ind.col] %>% as.matrix() #in case there is only one cell in the cluster
      temp <- softmax(as.vector(mat.cl))
      mat.cl.norm <- matrix(temp, nrow=nrow(mat.cl), ncol=ncol(mat.cl))
      # ## final re-scaling to 0 and 1 PER CLUSTER values
      # mat.cl.norm.rescale <- rescale(mat.cl.norm)
      # mat.norm[,ind.col] <- mat.cl.norm.rescale
      ## L2 normalization (i.e. unit vector length)
      mat.cl.norm.unit <- apply(mat.cl.norm, 2, FUN = function(z){z/sqrt(sum(z^2))})
      mat.norm[,ind.col] <- mat.cl.norm.unit
      # assay(sce.object) <- mat.norm
    }
    assay(sce.object) <- mat.norm
  }

  ## normalize by row mean PER CLUSTER
  if(norm.by=="mean"){
    mat.norm <- assay(sce.object)
    clusters <- unique(colData(sce.object)$cluster_membership)
    for(cl in clusters){
      ind.col <- colData(sce.object)$cluster_membership==cl
      mat.cl <- mat.norm[,ind.col] %>% as.matrix() #in case there is only one cell in the cluster
      mat.cl.norm <- sweep(mat.cl,1,rowMeans(mat.cl),"*")
      ## final re-scaling to 0 and 1 PER CLUSTER values
      mat.cl.norm.rescale <- mat.cl.norm/max(mat.cl.norm)
      mat.norm[,ind.col] <- mat.cl.norm.rescale
      # assay(sce.object) <- mat.norm
    }
    assay(sce.object) <- mat.norm
  }

  # ## normalize by row median PER CLUSTER
  # if(norm.by=="median"){
  #   mat.norm <- assay(sce.object)
  #   clusters <- unique(colData(sce.object)$cluster_membership)
  #   for(cl in clusters){
  #     ind.col <- colData(sce.object)$cluster_membership==cl
  #     mat.cl <- mat.norm[,ind.col] %>% as.matrix() #in case there is only one cell in the cluster
  #     mat.cl.norm <- sweep(mat.cl,1,rowMedians(mat.cl),"*")
  #     ## final re-scaling to 0 and 1 PER CLUSTER values
  #     mat.cl.norm.rescale <- mat.cl.norm/max(mat.cl.norm)
  #     mat.norm[,ind.col] <- mat.cl.norm.rescale
  #     # assay(sce.object) <- mat.norm
  #   }
  #   assay(sce.object) <- mat.norm
  # }

  ## do nothing
  if(norm.by=="none") sce.object <- sce.object

  return(sce.object)
}

