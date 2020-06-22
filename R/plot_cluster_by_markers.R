
#' "Barcoding" plot
#'
#' This function plots expression pattern of a cluster by marker genes, a.k.a. barcoding plot.
#'
#' @param sce.E1,sce.E2 Data objects of the \link[SingleCellExperiment]{SingleCellExperiment} data class. If only \code{sce.E1},
#' then the "barcode" is a self plot, i.e. both cluster (\code{cluster.name}) and marker genes are from the same experiment.
#' If both \code{sce.E1} and \code{sce.E2} are provided, then the "barcode" is a cross-experiment plot, i.e. marker genes are from
#' \code{sce.E1} (reference) andcluster (\code{cluster.name}) is from \code{sce.E2} (query).
#' @param cluster.name Name of the cluster to be plotted.
#' @param nsamp Number of randomly selected cells to plot for a large cluster. Default: \code{30}.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param use.common.markergenes Boolean variable indicating if to plot only common marker genes in a cross-experiment plot.
#' @param cellwidth,cellheight,main,filename,... Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @export

plot_cluster_by_markers <- function(sce.E1, sce.E2=NULL, cluster.name, nsamp=30,
                                    name.E1="E1", name.E2="E2", use.common.markergenes=TRUE,
                                    cellheight=10, cellwidth=5, main=NULL, filename=NA, ...){
  sce.ref <- sce.E1
  if(is.null(sce.E2)) sce.query <- sce.E1 else sce.query <- sce.E2

  ## check if the cluster is found in the sce.object
  if(!cluster.name %in% unique(colData(sce.query)$cluster_membership)){
    stop(paste(cluster.name, "is not found in the plotting data object. \n"))}

  ## reference marker genes
  markergenes <- rownames(sce.ref)[rowData(sce.ref)$marker_gene==1]
  ## cells of query cluster
  col.query <- colData(sce.query)$cluster_membership==cluster.name

  if(use.common.markergenes){
    ## matrix of common markergenes and query cells
    markergenes.common <- base::intersect(markergenes, rownames(sce.query))
    mat.query <- assay(sce.query[markergenes.common, col.query])
  } else{
    ## matrix of markergenes and query cells
    mat.query <- assay(sce.query[,col.query]) %>% as.data.frame() %>% rownames_to_column() %>%
      right_join(as.data.frame(markergenes, stringsAsFactors=FALSE), by=c("rowname"="markergenes")) %>%
      column_to_rownames() %>% as.matrix()
  }

  ## randomly select nsamp number of cells
  if(ncol(mat.query)>nsamp) mat.query <- mat.query[,sample(1:ncol(mat.query), nsamp, replace=FALSE)]

  ### self plot ###
  if(is.null(sce.E2)){
    ## main
    if(is.null(main)) main <- paste0(name.E1,".",cluster.name)
    ## indicator for markers
    if(!is.null(sce.query@metadata$cluster_marker_info)){
      mat.query <- mat.query[unique(sce.query@metadata$cluster_marker_info$markerGene),]
      temp <- rep(0, nrow(mat.query))
      markers.i <- sce.query@metadata$cluster_marker_info %>% filter(cluster==cluster.name) %>% pull(markerGene)
      temp[rownames(mat.query) %in% markers.i] <- 1
      ann <- data.frame("Marker"=as.factor(temp))
      rownames(ann) <- rownames(mat.query)
    } else {
      ann <- NA
    }
    ## plot
    pheatmap(mat.query,
             color = viridis::inferno(10), border_color = NA,
             cluster_rows = FALSE, cluster_cols = FALSE,
             cellheight = cellheight, cellwidth = cellwidth,
             show_rownames = TRUE, show_colnames = FALSE,
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
             color = viridis::inferno(10), border_color = NA,
             cluster_rows = FALSE, cluster_cols = FALSE,
             cellheight = cellheight, cellwidth = cellwidth,
             show_rownames = TRUE, show_colnames = FALSE,
             labels_col = "Cells", angle_col = "0",
             main=main,
             filename=filename,
             ...)
  }
}











