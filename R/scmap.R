
#' @export

######################
## match_by_scmap() ##
######################

match_by_scmap <- function(sce.query, sce.ref, imputation=FALSE,
                           filter.size=20, filter.fscore=NULL, #filter clusters as FRmatch
                           use.markers=FALSE, suppress_plot=TRUE){
  ## filtering small or low fscore clusters
  cat("Filtering small clusters: clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)) cat("Filtering low f-measure clusters in the reference experiment only. \n")
  sce.ref <- filter.cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)
  sce.query <- filter.cluster(sce.query, filter.size=filter.size, filter.fscore=NULL)

  ## imputation
  if(imputation){
    sce.ref <- impute_dropout(sce.ref)
    cat("Imputation done. \n")
  }

  ## extract info from sce.objects
  oo.query <- sce.query@metadata$cluster_order
  oo.ref <- sce.ref@metadata$cluster_order

  ## modify the sce.objects for scmap
  cat("Preparing data for scmap... \n")
  SummarizedExperiment::rowData(sce.query)$feature_symbol <- rownames(sce.query)
  SummarizedExperiment::colData(sce.query)$cell_type1 <- SummarizedExperiment::colData(sce.query)$cluster_membership
  SummarizedExperiment::rowData(sce.ref)$feature_symbol <- rownames(sce.ref)
  SummarizedExperiment::colData(sce.ref)$cell_type1 <- SummarizedExperiment::colData(sce.ref)$cluster_membership

  ## --------feature selection--------
  ## manually select features of NS-Forest marker genes
  if(use.markers){
    SummarizedExperiment::rowData(sce.ref)$scmap_features <- SummarizedExperiment::rowData(sce.ref)$NSF_markers
    cat("Customized markers are used. \n")
  }
  ## scmap selects 500 features
  else sce.ref <- scmap::selectFeatures(sce.ref, suppress_plot=suppress_plot)

  ## --------find median of features--------
  sce.ref <- scmap::indexCluster(sce.ref)

  ## scmap projection
  cat("Matching by scmapCluster... \n")
  scmapCluster_results <- scmap::scmapCluster(
    projection = sce.query,
    index_list = list(
      ref = metadata(sce.ref)$scmap_cluster_index
    )
  )

  ## matches
  prediction.match <- cbind(scmapCluster_results$scmap_cluster_labs[,'ref'],
                            SingleCellExperiment::colData(sce.query)$cell_type1)
  colnames(prediction.match) <- c("prediction", "celltype")

  celltypes <- unique(prediction.match[,"celltype"])
  match <- vector("list",length=length(celltypes))
  names(match) <- celltypes
  for(celltype in celltypes){
    ind <- prediction.match[,"celltype"]==celltype
    match[[celltype]] <- table(prediction.match[ind,1])
  }
  match <- match[oo.query]
  pctmat <- match2mat.scmap(match, oo.ref)

  ## return
  return(list("method"="scmapCluster", "matched.list"=match, "pctmat"=pctmat))
}

#######################
## match2mat.scmap() ##
#######################

match2mat.scmap <- function(match, oo.ref){
  df.match <- tibble::rownames_to_column(data.frame("freq"=unlist(match)))
  df.all <- merge(names(match), c(oo.ref,"unassigned"), all=TRUE) %>%
    dplyr::mutate(query.ref=paste0(x,".",y)) %>%
    dplyr::left_join(df.match, by=c("query.ref"="rowname")) %>%
    dplyr::mutate(freq=tidyr::replace_na(freq,0))
  mat <- matrix(df.all$freq, ncol=length(match), byrow=TRUE)
  rownames(mat) <- c(paste0("ref.",oo.ref), "unassigned") #row=ref
  colnames(mat) <- paste0("query.",names(match)) #col=query

  mat.pct <- sweep(mat,2,colSums(mat),"/")
  return(mat.pct)
}

##################
## plot_scmap() ##
##################

plot_scmap <- function(rst.scmap, type="matches", pct.cutoff=0.5,
                       reorder=TRUE, return.value=FALSE,
                       cellwidth=10, cellheight=10, main=NULL, ...){
  pctmat <- rst.scmap$pctmat
  pctmat.cutoff <- cutoff.pctmat(pctmat, pct.cutoff=pct.cutoff)
  if(reorder) pctmat.cutoff <- FRmatch:::reorder(pctmat.cutoff)
  if(type=="matches"){
    if(is.null(main)) main <- "Cell-to-cluster matches"
    pheatmap::pheatmap(pctmat.cutoff,
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")[c(3,3,7)]))(3),
                       breaks = seq(0,1,length.out=3),
                       legend_breaks=c(0,1), legend_labels=c("No match", "Match"),
                       cluster_rows=F, cluster_cols=F,
                       gaps_row=nrow(pctmat.cutoff)-1,
                       cellwidth=cellwidth, cellheight=cellheight,
                       main=main,
                       ...)
    if(return.value) return(pctmat.cutoff)
  }
  if(type=="pctmat"){
    if(reorder) pctmat <- pctmat[,colnames(pctmat.cutoff)]
    if(is.null(main)) main <- "Percent of cells matched"
    pheatmap::pheatmap(pctmat,
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101),
                       breaks = seq(0,1,length.out=101),
                       cluster_rows=F, cluster_cols=F,
                       cellwidth=cellwidth, cellheight=cellheight,
                       gaps_row=nrow(pctmat)-1,
                       main=main,
                       ...)
    if(return.value) return(pctmat)
  }
}

cutoff.pctmat <- function(pctmat, pct.cutoff){
  pctmat.cutoff <- matrix(as.numeric(pctmat>=pct.cutoff), nrow=nrow(pctmat))
  pctmat.cutoff[nrow(pctmat.cutoff),] <- as.numeric(colSums(pctmat.cutoff[-nrow(pctmat.cutoff),])==0)
  rownames(pctmat.cutoff) <- rownames(pctmat)
  colnames(pctmat.cutoff) <- colnames(pctmat)
  return(pctmat.cutoff)
}

############################
## plot_bilateral_scmap() ##
############################

plot_bilateral_scmap <- function(rst.scmap.E1toE2, rst.scmap.E2toE1, prefix.E1="E1", prefix.E2="E2",
                                 pct.cutoff=0.5,
                                 reorder=TRUE, return.value=FALSE,
                                 cellwidth=10, cellheight=10, main=NULL, ...){
  ## get binary matrices for plotting
  pctmat.cutoff.E1toE2 <- cutoff.pctmat(rst.scmap.E1toE2$pctmat, pct.cutoff=pct.cutoff)
  pctmat.cutoff.E2toE1 <- cutoff.pctmat(rst.scmap.E2toE1$pctmat, pct.cutoff=pct.cutoff)

  ## combine two matrices to one bilateral matrix
  mat1 <- pctmat.cutoff.E1toE2[-nrow(pctmat.cutoff.E1toE2),] #use E1toE2 as the framework for final plot
  mat2 <- t(pctmat.cutoff.E2toE1[-nrow(pctmat.cutoff.E2toE1),]) #so transpose E2toE1
  mat.bi <- mat1+mat2
  ## unassigned row
  mat.bi <- rbind(mat.bi, 2*as.numeric(colSums(mat.bi)==0))
  ## rename colnames and rownames
  rownames(mat.bi) <- gsub("ref",prefix.E2,rownames(pctmat.cutoff.E1toE2))
  colnames(mat.bi) <- gsub("query",prefix.E1,colnames(pctmat.cutoff.E1toE2))

  ## plot
  if(is.null(main)) main <- "Bilatreral matches (cell-to-cluster)"
  if(reorder) mat.bi <- FRmatch:::reorder(mat.bi)
  pheatmap::pheatmap(mat.bi,
                     color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")[c(1,1,3,7)]))(4),
                     breaks = seq(0,2,length.out=4),
                     legend_breaks=0:2, legend_labels=c("No match", "One-way match", "Two-way Match"),
                     cluster_rows=F, cluster_cols=F,
                     gaps_row=nrow(mat.bi)-1,
                     cellwidth=cellwidth, cellheight=cellheight,
                     main=main,
                     ...)
  if(return.value) return(mat.bi)
}
