
#' Plot a cluster in the dimension of reference markers (barcoding plot)
#'
#' This is an auxilirary function that plots the expression levels of a (query) cluster by reference marker genes, a.k.a. barcoding plot.
#'
#' @param sce.E1,sce.E2 Data objects of the \link[SingleCellExperiment]{SingleCellExperiment} data class. If only \code{sce.E1} is provided,
#' it serves as both the query and reference dataset, which means that \code{cluster.name} should be one of the clusters in \code{sce.E1}
#' and it should also have marker genes. If both \code{sce.E1} and \code{sce.E2} are provided, then \code{sce.E1} will serve as the
#' reference dataset and \code{sce.E2} as the query dataset, which means that \code{cluster.name} should be one of the clusters in
#' \code{sce.E2} and the cluster will be plotted using marker genes from \code{sce.E1}.
#' @param cluster.name Name of the cluster to be plotted.
#' @param nsamp Number of randomly selected cells to plot for a large cluster. Default: \code{30}.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param cellwidth,cellheight,main,filename,... Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @export

plot_cluster_by_markers <- function(sce.E1, sce.E2=NULL, cluster.name, nsamp=30,
                                    name.E1="E1", name.E2="E2",
                                    cellheight=10, cellwidth=5, main=NULL, filename=NA, ...){
  sce.ref <- sce.E1
  if(is.null(sce.E2)) sce.query <- sce.E1 else sce.query <- sce.E2

  ## check if the cluster is found in the sce.object
  if(!cluster.name %in% unique(colData(sce.query)$cluster_membership)){
    stop(paste(cluster.name,
               "is not found in the SingleCellExperiment (sce) data object. \n"))}

  ## reference marker genes
  markergenes <- sce.ref@metadata$cluster_marker_info$markerGene
  markergenes.common <- base::intersect(markergenes, rownames(sce.query))

  ## query cluster
  col.query <- colData(sce.query)$cluster_membership==cluster.name
  mat.query <- counts(sce.query[markergenes.common, col.query])

  ## randomly select nsamp number of cells
  if(ncol(mat.query)>nsamp) mat.query <- mat.query[,sample(1:ncol(mat.query), nsamp, replace=FALSE)]

  ### self plot ###
  if(is.null(sce.E2)){
    ## main
    if(is.null(main)) main <- paste0(name.E1,".",cluster.name)
    ## indicator for markers
    temp <- rep(0, nrow(mat.query))
    markers.i <- sce.query@metadata$cluster_marker_info %>% filter(cluster==cluster.name) %>% pull(markerGene)
    temp[rownames(mat.query) %in% markers.i] <- 1
    ann <- data.frame("Marker"=as.factor(temp))
    rownames(ann) <- rownames(mat.query)
    ## plot
    pheatmap(mat.query,
             cluster_rows = F, cluster_cols = F,
             cellheight = cellheight, cellwidth = cellwidth,
             show_colnames = TRUE,
             annotation_row = ann,
             labels_col = "Cells", angle_col = "0",
             main=main,
             filename=filename,
             ...)
  }

  ### cross-experiment plot ###
  if(!is.null(sce.E2)){
    ## main
    if(is.null(main)) main <- paste0(name.E2,".",cluster.name)
    ## plot
    pheatmap(mat.query,
             cluster_rows = F, cluster_cols = F,
             cellheight = cellheight, cellwidth = cellwidth,
             show_colnames = TRUE,
             labels_col = "Cells", angle_col = "0",
             main=main,
             filename=filename,
             ...)
  }
}











