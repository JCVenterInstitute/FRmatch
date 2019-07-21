


#' @export

#######################
## match_by_seurat() ##
#######################

match_by_seurat <- function(sce.query, sce.ref, imputation=FALSE,
                            filter.size=20, filter.fscore=NULL, #filter clusters as FRmatch
                            use.markers=FALSE){

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

  ## make data into seurat objects in a list
  cat("Preparing data for Seurat... \n")
  mat <- merge(SummarizedExperiment::assay(sce.query), SummarizedExperiment::assay(sce.ref), by="row.names")
  rownames(mat) <- mat$Row.names
  mat <- as.matrix(mat[,-1])
  data <- as(mat, "dgCMatrix")

  metadata <- data.frame("subsample"=rep(c("query","ref"), c(ncol(sce.query),ncol(sce.ref))),
                         "celltype"=c(as.character(SingleCellExperiment::colData(sce.query)$cluster_membership),
                                      as.character(SingleCellExperiment::colData(sce.ref)$cluster_membership)),
                         stringsAsFactors=FALSE)
  rownames(metadata) <- colnames(mat)

  seurat.object <- Seurat::CreateSeuratObject(counts=data, meta.data=metadata)
  seurat.object.list <- Seurat::SplitObject(object=seurat.object, split.by="subsample")

  ## select variable features
  for (i in 1:length(x=seurat.object.list)) {
    seurat.object.list[[i]] <- Seurat::FindVariableFeatures(object=seurat.object.list[[i]],
                                                            selection.method="vst", nfeatures=2000, verbose=FALSE)
  }

  ### matching ###
  cat("Matching by Seurat... \n")
  object.ref <- seurat.object.list[[c("ref")]]
  object.query <- seurat.object.list[[c("query")]]
  ## PCA used by defualt in the following step, which is recommended for scRNA-seq
  anchors <- Seurat::FindTransferAnchors(reference=object.ref, query=object.query, dims=1:30)
  predictions <- Seurat::TransferData(anchorset=anchors, refdata=object.ref$celltype, dims=1:30)

  ## calculate Seurat scores averaged by cluster
  df.predictions <- predictions %>% tibble::rownames_to_column()
  df.predictions$rowname <- gsub(".x","",df.predictions$rowname)
  df.membership <- SingleCellExperiment::colData(sce.query) %>% as.data.frame() %>% tibble::rownames_to_column()
  df.scores <- df.predictions %>% dplyr::left_join(df.membership) %>% dplyr::group_by(cluster_membership) %>%
    dplyr::select(-c(rowname, predicted.id, prediction.score.max)) %>%
    dplyr::summarise_all(mean) # average scores by cluster
  mymat.scores <- df.scores %>% dplyr::select(dplyr::starts_with("prediction")) %>% as.matrix() %>% t()
  colnames(mymat.scores) <- paste0("query.", df.scores$cluster_membership)
  scores <- mymat.scores[paste0("prediction.score.",oo.ref),paste0("query.",oo.query)]

  object.query <- Seurat::AddMetaData(object=object.query, metadata=predictions)
  prediction.match <- cbind(object.query$predicted.id, object.query$celltype)
  colnames(prediction.match) <- c("prediction", "celltype")

  celltypes <- unique(prediction.match[,"celltype"])
  match <- vector("list",length=length(celltypes))
  names(match) <- celltypes
  for(celltype in celltypes){
    ind <- prediction.match[,"celltype"]==celltype
    match[[celltype]] <- table(prediction.match[ind,1])
  }
  match <- match[oo.query]

  cluster.sizes.query <- table(SingleCellExperiment::colData(sce.query)$cluster_membership)[names(match)]
  pctmat <- match2mat.seurat(match, oo.ref, cluster.sizes.query)

  ## return
  return(list("pctmat"=pctmat, "prediction.scores"=scores))
}

########################
## match2mat.seurat() ##
########################

match2mat.seurat <- function(match, oo.ref, cluster.sizes.query){
  df.match <- tibble::rownames_to_column(data.frame("freq"=unlist(match)))
  df.all <- merge(names(match), oo.ref, all=TRUE) %>%
    dplyr::mutate(query.ref=paste0(x,".",y)) %>%
    dplyr::left_join(df.match, by=c("query.ref"="rowname")) %>%
    dplyr::mutate(freq=tidyr::replace_na(freq,0))
  mat <- matrix(df.all$freq, nrow=length(oo.ref), byrow=TRUE)

  # cluster.sizes.query <- table(SingleCellExperiment::colData(sce.query)$cluster_membership)[names(match)]
  unassigned <- cluster.sizes.query - colSums(mat)
  mat <- rbind(mat,unassigned)

  rownames(mat) <- c(paste0("ref.",oo.ref), "unassigned") #row=ref
  colnames(mat) <- paste0("query.",names(match)) #col=query

  mat.pct <- sweep(mat,2,cluster.sizes.query,"/")
  return(mat.pct)
}

###################
## plot_seurat() ##
###################

plot_seurat <- function(rst.seurat, type="matches", pct.cutoff=0.5,
                        reorder=TRUE, return.value=FALSE,
                        cellwidth=10, cellheight=10, main=NULL, ...){
  pctmat <- rst.seurat$pctmat
  pctmat.cutoff <- FRmatch:::cutoff.pctmat(pctmat, pct.cutoff=pct.cutoff)
  if(reorder) pctmat.cutoff <- FRmatch:::reorder(pctmat.cutoff)
  if(type=="matches"){
    if(is.null(main)) main <- "Cell-to-cell matches"
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
  if(type=="prediction.scores"){
    prediction.scores <- rst.seurat$prediction.scores
    if(reorder) prediction.scores <- prediction.scores[,colnames(pctmat.cutoff)]
    if(is.null(main)) main <- "Prediction score"
    pheatmap::pheatmap(prediction.scores,
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(101),
                       breaks = seq(0,1,length.out=101),
                       cluster_rows=F, cluster_cols=F,
                       cellwidth=cellwidth, cellheight=cellheight,
                       main=main,
                       ...)
    if(return.value) return(prediction.scores)
  }
}

#############################
## plot_bilateral_seurat() ##
#############################

plot_bilateral_seurat <- function(rst.seurat.E1toE2, rst.seurat.E2toE1, prefix.E1="E1", prefix.E2="E2",
                                 pct.cutoff=0.5,
                                 reorder=TRUE, return.value=FALSE,
                                 cellwidth=10, cellheight=10, main=NULL, ...){
  ## get binary matrices for plotting
  pctmat.cutoff.E1toE2 <- cutoff.pctmat(rst.seurat.E1toE2$pctmat, pct.cutoff=pct.cutoff)
  pctmat.cutoff.E2toE1 <- cutoff.pctmat(rst.seurat.E2toE1$pctmat, pct.cutoff=pct.cutoff)

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
  if(is.null(main)) main <- "Bilatreral matches (cell-to-cell)"
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
