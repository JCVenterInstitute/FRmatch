## reqires
# pheatmap::pheatmap
# RColorBrewer::brewer.pal


#' Check per-cluster dropout rate for cluster marker genes
#'
#' A function that calculates and plots the dropout rate for a \code{SingleCellExperiment} object
#' customized with cluster and NS-Forest marker gene information.
#'
#' @param sce.object A \code{SingleCellExperiment} object filled with necessary information for \code{FRmatch}.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param return.value Logical variable indicating if to return the dropout rate values to R console. Default: \code{TRUE}.
#' @param plot.dropout Logical variable indicating if to return a plot of dropout rates. Default: \code{FALSE}.
#' @param cellwidth,cellheight,... Additional plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return A numeric matrix of dropout rates per cluster per marker gene, and optionally, a plot.
#'
#' @author Yun Zhang, \email{zhangy@jcvi.org}; Brian Aevermann, \email{baeverma@jcvi.org}; Richard Scheuermann, \email{RScheuermann@jcvi.org}.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' check_dropout(sce.example, plot.dropout=TRUE, main="Dropout rate plot")
#' }
#'
#' @importFrom dplyr %>%
#' @export

check_dropout <- function(sce.object, return.value=TRUE, plot.dropout=FALSE, cellwidth=15, cellheight=10, ...){
  ## data
  dat <- logcounts(sce.object)
  fscores <- sce.object@metadata$fscores

  ## cluster info
  cluster_order <- sce.object@metadata$cluster_order
  cluster_marker_info <- metadata(sce.object)$cluster_marker_info %>% dplyr::arrange(match(cluster, cluster_order))
  cluster_membership <- colData(sce.object)$cluster_membership

  ## pct of zeros per marker gene per cluster
  cluster_marker_mat <- dat[cluster_marker_info$markerGene,]
  dat.lst <- lapply(split(as.data.frame(t(cluster_marker_mat)), cluster_membership),t)
  dat.lst <- dat.lst[cluster_order]
  zero.pct <- sapply(dat.lst, function(z) rowSums(z==0)/ncol(z))

  ## plot
  if(plot.dropout){
    gaps <- cumsum(table(cluster_marker_info$cluster)[cluster_order])
    if(!is.null(fscores)){
      ann <- data.frame("f.score"=fscores$`f-measure`)
      rownames(ann) <- fscores$clusterName
      pheatmap::pheatmap(zero.pct,
                         color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Blues")))(101),
                         breaks=seq(0,1,length.out=101),
                         cluster_rows = F, cluster_cols = F,
                         gaps_row = gaps[-length(gaps)], gaps_col = 1:(ncol(zero.pct)-1),
                         annotation_col=ann,
                         cellwidth=cellwidth, cellheight=cellheight, ...)
    }
    else
    pheatmap::pheatmap(zero.pct,
                       color=colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "Blues")))(101),
                       breaks=seq(0,1,length.out=101),
                       cluster_rows = F, cluster_cols = F,
                       gaps_row = gaps[-length(gaps)], gaps_col = 1:(ncol(zero.pct)-1),
                       cellwidth=cellwidth, cellheight=cellheight, ...)
  }

  ## output
  if(return.value) return(zero.pct)
}
