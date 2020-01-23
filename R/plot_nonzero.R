
#' Non-zero expression plot
#'
#' A function that calculates and plots "\% expressed per marker gene per cluster" for the \code{FRmatch} input data object.
#' \% = number of cells that express the marker gene in the cluster / cluster size.
#'
#' @param sce.object A \code{FRmatch} input data object. See example in \code{\link[FRmatch]{sce.example}}.
#' @param return.plot Logical variable indicating if to return the plot. Default: \code{TRUE}.
#' @param return.value Logical variable indicating if to return the plotted values. Default: \code{FALSE}.
#' @param cellwidth,cellheight,main,... Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return if \code{return.value = TRUE}, a matrix of plotted values.
#'
#' @examples
#' \dontrun{
#' data("sce.example")
#' plot_nonzero(sce.example)
#' }
#'
#' @export

plot_nonzero <- function(sce.object, return.plot=TRUE, return.value=FALSE,
                         cellwidth=15, cellheight=10, main=NULL, ...){
  ## data
  dat <- counts(sce.object)
  fscores <- sce.object@metadata$fscores

  ## cluster info
  cluster_order <- sce.object@metadata$cluster_order
  cluster_marker_info <- sce.object@metadata$cluster_marker_info %>% arrange(match(cluster, cluster_order))
  cluster_membership <- colData(sce.object)$cluster_membership

  ## pct of zeros per marker gene per cluster
  cluster_marker_mat <- dat[cluster_marker_info$markerGene,]
  dat.lst <- lapply(split(as.data.frame(t(cluster_marker_mat)), cluster_membership),t)
  dat.lst <- dat.lst[cluster_order]
  zero.pct <- sapply(dat.lst, function(z) rowSums(z==0)/ncol(z))
  nonzero.pct <- 1-zero.pct

  ## plot
  if(return.plot){
    if(is.null(main)) main <- "% expressed per marker gene per cluster"
    gaps <- cumsum(table(cluster_marker_info$cluster)[cluster_order])
    if(!is.null(fscores)){
      ann <- data.frame("f.score"=fscores$`f-measure`)
      rownames(ann) <- fscores$clusterName
      pheatmap(nonzero.pct,
               color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101),
               breaks=seq(0,1,length.out=101),
               cluster_rows = F, cluster_cols = F,
               gaps_row = gaps[-length(gaps)], gaps_col = 1:(ncol(zero.pct)-1),
               annotation_col=ann,
               cellwidth=cellwidth, cellheight=cellheight,
               main=main, ...)
    }
    else
      pheatmap(nonzero.pct,
               color=colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101),
               breaks=seq(0,1,length.out=101),
               cluster_rows = F, cluster_cols = F,
               gaps_row = gaps[-length(gaps)], gaps_col = 1:(ncol(zero.pct)-1),
               cellwidth=cellwidth, cellheight=cellheight,
               main=main, ...)
  }

  ## output
  if(return.value) return(nonzero.pct)
}
