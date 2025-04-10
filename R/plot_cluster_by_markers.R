
#' "Barcode" plot
#'
#' This function plots expression pattern of a cluster by marker genes, a.k.a. barcode plot.
#'
#' @param sce.E1,sce.E2 Data objects of the \link[SingleCellExperiment]{SingleCellExperiment} data class. If only \code{sce.E1},
#' then the "barcode" is a self plot, i.e. both cluster (\code{cluster.name}) and marker genes are from the same experiment.
#' If both \code{sce.E1} and \code{sce.E2} are provided, then the "barcode" is a cross-experiment plot, i.e. marker genes are from
#' \code{sce.E1} (reference) and cluster (\code{cluster.name}) is from \code{sce.E2} (query).
#' @param cluster.name Name of the cluster to be plotted.
#' @param nsamp Number of randomly selected cells to plot for a large cluster. Default: \code{30}.
#' @param name.self,name.cross Prefix name of experiment for self and cross-experiment plots. Default: \code{"E1."} and \code{"E2."}, respectively.
#' @param use.common.markergenes Boolean variable indicating if to plot only common marker genes in a cross-experiment plot.
#' @param scale.colorbar Boolean variable indicating if to scale the color bar to [0,1] for normalized gene expression values. Default: \code{"FALSE"}.
#' @param cellwidth,cellheight,main,filename,... Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @import viridis
#'
#' @export

plot_cluster_by_markers <- function(sce.E1, sce.E2=NULL, cluster.name, nsamp=30,
                                    name.self="E1.", name.cross="E2.", use.common.markergenes=TRUE,
                                    scale.colorbar=FALSE, cellheight=10, cellwidth=5, main=NULL, filename=NA, ...){
  sce.ref <- sce.E1
  if(is.null(sce.E2)) sce.query <- sce.E1 else sce.query <- sce.E2

  ## check if the cluster is found in the sce.object
  if(!cluster.name %in% unique(colData(sce.query)$cluster_membership)){
    stop(paste(cluster.name, "is not found in the plotting data object. \n"))}

  ## REORDER clusters according to the given order if available
  if(!is.null(sce.ref@metadata$cluster_order)){
    sce.ref@metadata$cluster_marker_info %<>% arrange(match(clusterName, sce.ref@metadata$cluster_order))}

  ## reference marker genes
  markergenes <- unique(sce.ref@metadata$cluster_marker_info$markerGene) #marker genes in ORDER!!!
  if(is.null(markergenes)) markergenes <- rownames(sce.ref)[rowData(sce.ref)$marker_gene==1] #if metadat is not available
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
    mat.query <- mat.query[markergenes,]
  }

  ## randomly select nsamp number of cells
  if(ncol(mat.query)>nsamp) mat.query <- mat.query[,sample(1:ncol(mat.query), nsamp, replace=FALSE)]

  ## if to scale colorbar to [0,1] for normalized expr
  if(scale.colorbar){
    breaks <- seq(0, 1, length.out = 11)
  } else breaks <- NA

  ### self plot ###
  if(is.null(sce.E2)){
    ## main
    if(is.null(main)) main <- paste0(name.self, cluster.name)
    ## indicator for markers
    if(!is.null(sce.query@metadata$cluster_marker_info)){
      # mat.query <- mat.query[unique(sce.query@metadata$cluster_marker_info$markerGene),]
      temp <- rep(0, nrow(mat.query))
      markers.i <- sce.query@metadata$cluster_marker_info %>% filter(clusterName==cluster.name) %>% pull(markerGene)
      temp[rownames(mat.query) %in% markers.i] <- 1
      ann <- data.frame("Marker"=as.factor(temp))
      rownames(ann) <- rownames(mat.query)
    } else {
      ann <- NA
    }
    ann_colors = list(
      Marker = c("0"="lightgrey", "1"="#00BFC4")
    )
    if(!is.null(sce.ref@metadata$f_score)){
      fscore <- sce.ref@metadata$f_score %>% filter(clusterName==cluster.name) %>% pull(score)
    }
    else fscore <- NA
    ## plot
    pheatmap(mat.query,
             color = inferno(10), breaks = breaks, border_color = NA,
             cluster_rows = FALSE, cluster_cols = FALSE,
             cellheight = cellheight, cellwidth = cellwidth,
             show_rownames = TRUE, show_colnames = TRUE, angle_col = "0",
             annotation_row = ann, annotation_colors=ann_colors,
             labels_col = paste("F-beta:", round(fscore,3)), fontsize_col = 15,
             main = main, filename = filename,
             ...)
  }

  ### cross-experiment plot ###
  if(!is.null(sce.E2)){
    ## main
    if(is.null(main)) main <- paste0(name.cross,cluster.name)
    ## plot
    pheatmap(mat.query,
             color = inferno(10), breaks = breaks, border_color = NA,
             cluster_rows = FALSE, cluster_cols = FALSE,
             cellheight = cellheight, cellwidth = cellwidth,
             show_rownames = TRUE, show_colnames = FALSE,
             labels_col = "Cells", angle_col = "0",
             main=main, filename=filename,
             ...)
  }
}











