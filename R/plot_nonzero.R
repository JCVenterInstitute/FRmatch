
#' Plot per-cluster percentage of non-zero expression for markers
#'
#' A function that calculates and plots percentage of non-zero expression for each marker gene per cluster
#' for a \code{SingleCellExperiment} object customized with cluster membership and marker gene information.
#'
#' @param sce.object A \code{SingleCellExperiment} object customized with necessary information for \code{FRmatch}.
#' See details in \code{\link[FRmatch]{sce.example}}.
#' @param return.plot Logical variable indicating if to return a plot. Default: \code{TRUE}.
#' @param return.value Logical variable indicating if to return the values correponding to the plot. Default: \code{FALSE}.
#' @param cellwidth,cellheight,main,... Additional plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return Optionally, a numeric matrix corresponding to the values on the plot.
#'
#' @seealso The \link[SingleCellExperiment]{SingleCellExperiment} class.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' plot_nonzero(sce.example)
#' }
#'
#' @importFrom dplyr %>%
#' @export

plot_nonzero <- function(sce.object, return.plot=TRUE, return.value=FALSE,
                         cellwidth=15, cellheight=10, main=NULL, ...){
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
  nonzero.pct <- 1-zero.pct

  ## plot
  if(return.plot){
    if(is.null(main)) main <- "% expressed per marker per cluster"
    gaps <- cumsum(table(cluster_marker_info$cluster)[cluster_order])
    if(!is.null(fscores)){
      ann <- data.frame("f.score"=fscores$`f-measure`)
      rownames(ann) <- fscores$clusterName
      pheatmap::pheatmap(nonzero.pct,
                         color=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Blues"))(101),
                         breaks=seq(0,1,length.out=101),
                         cluster_rows = F, cluster_cols = F,
                         gaps_row = gaps[-length(gaps)], gaps_col = 1:(ncol(zero.pct)-1),
                         annotation_col=ann,
                         cellwidth=cellwidth, cellheight=cellheight,
                         main=main, ...)
    }
    else
    pheatmap::pheatmap(nonzero.pct,
                       color=colorRampPalette(RColorBrewer::brewer.pal(n = 7, name = "Blues"))(101),
                       breaks=seq(0,1,length.out=101),
                       cluster_rows = F, cluster_cols = F,
                       gaps_row = gaps[-length(gaps)], gaps_col = 1:(ncol(zero.pct)-1),
                       cellwidth=cellwidth, cellheight=cellheight,
                       main=main, ...)
  }

  ## output
  if(return.value) return(nonzero.pct)
}
