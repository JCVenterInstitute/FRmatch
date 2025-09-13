
#' Plot minimum spanning tree (MST)
#'
#' This function is a wrapper function for plotting MST of two interested clusters.
#'
#' @param sce.query,sce.ref Query and reference data objects.
#' @param query.cluster,ref.cluster Query and reference cluster names to plot.
#' @param nsamp Number of randomly selected cells to plot for large cluster. Default: 30.
#' @param ... Additional arguments passed to \code{\link[FRmatch]{FRtest}}.
#'
#' @return MST plot and FR-test result in console.
#'
#' @export

plot_MST <- function(sce.query, sce.ref, query.cluster, ref.cluster, nsamp=30, ...){
  ind.query <- sce.query@colData$cluster_membership==query.cluster
  ind.query.sub <- sample(1:sum(ind.query), min(nsamp,sum(ind.query)))
  samp1 <- assay(sce.query)[,ind.query][,ind.query.sub]

  ind.ref <- sce.ref@colData$cluster_membership==ref.cluster
  ind.ref.sub <- sample(1:sum(ind.ref), min(nsamp,sum(ind.ref)))
  samp2 <- assay(sce.ref)[,ind.ref][,ind.ref.sub]

  FRtest(samp1, samp2, plot.MST=T, label.names=c(query.cluster,ref.cluster))
}
