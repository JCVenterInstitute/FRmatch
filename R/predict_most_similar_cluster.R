
#' Prediction function for FR-Match results
#'
#' This function takes in the \code{\link[FRmatch]{FRmatch}} output and predicts the most similar reference cluster for each query cluster.
#'
#' @param rst.FRmatch The \code{\link[FRmatch]{FRmatch}} output.
#' @param p.adj.method P-value adjustment method for multiple hypothesis testing correction. Default: \code{"BY"}.
#' For more options, please see \code{\link[stats]{p.adjust.methods}}.
#'
#' @return A data frame of the most similar reference cluster for each query cluster,
#' which is the "match" with the highest adjusted p-value.
#'
#' @export

predict_most_similar_cluster <- function(rst.FRmatch, p.adj.method="BY"){
  pmat <- rst.FRmatch$pmat
  pmat.adj <- padj.FRmatch(pmat, p.adj.method=p.adj.method)
  df.most.similar.cluster <- data.frame("query_cluster"=colnames(pmat),
                                        "most_similar_ref_cluster"=rownames(pmat)[apply(pmat, 2, which.max)], #use the original pmat, independent of p.adj.method
                                        "padj"=colMaxs(pmat.adj)) %>% arrange(desc(padj)) #return padj here
  return(df.most.similar.cluster)
}
