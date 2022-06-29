
# #' Functions for other matching methods: scmap, Seurat
# #'
# #' @return A list of:
# #' \item{method}{Call of matching method used.}
# #' \item{pctmat}{A matrix of % matched cells. Rows are reference clusters, and columns are query clusters.}
# #' \item{prediction.scores}{A matrix of Seurat prediction scores. Rows are reference clusters, and columns are query clusters.}
# #' \item{all.results}{A data frame with columns \code{query.cell, query.cluster, predicted.ref.cluster}}.
# #'
# #' @export match_by_scmap,match_by_seurat #how to export multiple functions???

######################################################################################################################
### scmap ############################################################################################################
######################################################################################################################

######################
## match_by_scmap() ##
######################

match_by_scmap <- function(sce.query, sce.ref, #imputation=FALSE,
                           filter.size=10, filter.fscore=NULL, #filter clusters
                           use.markers=FALSE, suppress_plot=TRUE){
  ## filtering small or low fscore clusters
  cat("Filtering small clusters: clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)) cat("Filtering low F-beta score clusters in the reference experiment only. \n")
  sce.ref <- filter_cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)
  sce.query <- filter_cluster(sce.query, filter.size=filter.size, filter.fscore=NULL)

  ## imputation
  # if(imputation){
  #   sce.ref <- impute_dropout(sce.ref)
  #   cat("Imputation done. \n")
  # }

  ## extract info from sce.objects
  oo.query <- sce.query@metadata$cluster_order
  oo.ref <- sce.ref@metadata$cluster_order

  ## modify the sce.objects for scmap
  cat("Preparing data for scmap... \n")
  rowData(sce.query)$feature_symbol <- rownames(sce.query)
  colData(sce.query)$cell_type1 <- colData(sce.query)$cluster_membership
  rowData(sce.ref)$feature_symbol <- rownames(sce.ref)
  colData(sce.ref)$cell_type1 <- colData(sce.ref)$cluster_membership

  ## --------feature selection--------
  ## manually select features of NS-Forest marker genes
  if(use.markers){
    rowData(sce.ref)$scmap_features <- rowData(sce.ref)$NSF_markers
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
      ref = sce.ref@metadata$scmap_cluster_index
    )
  )

  ## prediction.celltype
  prediction.celltype <- cbind(scmapCluster_results$scmap_cluster_labs[,'ref'],
                               colData(sce.query)$cell_type1)
  colnames(prediction.celltype) <- c("prediction", "celltype")
  rownames(prediction.celltype) <- colnames(sce.query)

  ## match.lst
  celltypes <- unique(prediction.celltype[,"celltype"])
  match.lst <- vector("list",length=length(celltypes))
  names(match.lst) <- celltypes
  for(celltype in celltypes){
    ind <- prediction.celltype[,"celltype"]==celltype
    match.lst[[celltype]] <- table(prediction.celltype[ind,1])
  }
  match.lst <- match.lst[oo.query]
  pctmat <- match2mat.scmap(match.lst, oo.ref)

  ## all results
  all.results <- data.frame("query.cell"=rownames(prediction.celltype), "query.cluster"=prediction.celltype[,"celltype"],
                            "predicted.ref.cluster"=prediction.celltype[,"prediction"])

  ## return
  return(list("method"="scmapCluster", "pctmat"=pctmat, "all.results"=all.results))
}

##-------------------##
## match2mat.scmap() ##
##-------------------##

match2mat.scmap <- function(match.lst, oo.ref){
  df.match <- rownames_to_column(data.frame("freq"=unlist(match.lst)))
  df.all <- merge(names(match.lst), c(oo.ref,"unassigned"), all=TRUE) %>%
    mutate(query.ref=paste0(x,".",y)) %>%
    left_join(df.match, by=c("query.ref"="rowname")) %>%
    mutate(freq=replace_na(freq,0))
  mat <- matrix(df.all$freq, ncol=length(match.lst), byrow=TRUE)
  rownames(mat) <- c(paste0("ref.",oo.ref), "unassigned") #row=ref
  colnames(mat) <- paste0("query.",names(match.lst)) #col=query

  mat.pct <- sweep(mat,2,colSums(mat),"/")
  return(mat.pct)
}

##################
## plot_scmap() ##
##################

plot_scmap <- function(rst.scmap, type="matches", pct.cutoff=0.3,
                       reorder=TRUE, return.value=FALSE,
                       cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){
  pctmat <- rst.scmap$pctmat
  pctmat.cutoff <- cutoff.pctmat(pctmat, pct.cutoff=pct.cutoff)
  if(reorder){
    pctmat.cutoff <- reorder(pctmat.cutoff)
    pctmat <- pctmat[,colnames(pctmat.cutoff)]
  }
  if(type=="matches"){
    if(is.null(main)) main <- "scmap"
    ## heatmap
    pheatmap(pctmat.cutoff,
             color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(3,3,7)]))(3),
             breaks=seq(0,1,length.out=3),
             legend_breaks=c(0,1), legend_labels=c("No match", "Match"),
             cluster_rows=F, cluster_cols=F,
             gaps_row=nrow(pctmat.cutoff)-1,
             cellwidth=cellwidth, cellheight=cellheight,
             main=main, filename=filename,
             ...)
    cat("pct.cutoff = ", pct.cutoff, "\n")
    ## output
    if(return.value) return(pctmat.cutoff)
  }
  if(type=="pctmat"){
    # if(is.null(main)) main <- "Distribution of % matched cells"
    pctmat <- pctmat[-nrow(pctmat),]
    df <- tibble(pct=as.vector(pctmat),
                 query_cluster = rep(colnames(pctmat), each=nrow(pctmat))) %>%
      mutate(query_cluster = fct_relevel(query_cluster, colnames(pctmat)))
    g <- ggplot(df, aes(x=query_cluster, y=pct)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 270, hjust = 0)) +
      geom_hline(linetype = "dashed", yintercept = pct.cutoff, color = "red") +
      scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      xlab("Query cluster") + ylab("% matched cells")
    ## save plot or plot on device
    if(!is.na(filename)) ggsave(filename, g, width=ncol(pctmat)*.2, height=5)
    else plot(g)
    cat("pct.cutoff = ", pct.cutoff, "\n")
    ## output
    if(return.value) return(pctmat)
  }
}

##-----------------##
## cutoff.pctmat() ##
##-----------------##
cutoff.pctmat <- function(pctmat, pct.cutoff){
  pctmat.cutoff <- matrix(as.numeric(pctmat>=pct.cutoff), nrow=nrow(pctmat))
  pctmat.cutoff[nrow(pctmat.cutoff),] <- as.numeric(colSums(pctmat.cutoff[-nrow(pctmat.cutoff),])==0)
  rownames(pctmat.cutoff) <- rownames(pctmat)
  colnames(pctmat.cutoff) <- colnames(pctmat)
  return(pctmat.cutoff)
}

#####################
## plot_bi_scmap() ##
#####################

plot_bi_scmap <- function(rst.scmap.E1toE2, rst.scmap.E2toE1, name.E1="E1", name.E2="E2",
                          pct.cutoff=0.3,
                          reorder=TRUE, return.value=FALSE,
                          cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){
  ## get binary matrices for plotting
  pctmat.cutoff.E1toE2 <- cutoff.pctmat(rst.scmap.E1toE2$pctmat, pct.cutoff=pct.cutoff)
  pctmat.cutoff.E2toE1 <- cutoff.pctmat(rst.scmap.E2toE1$pctmat, pct.cutoff=pct.cutoff)

  ## combine two matrices to one two-way matrix
  mat1 <- pctmat.cutoff.E1toE2[-nrow(pctmat.cutoff.E1toE2),] #use E1toE2 as the framework for final plot
  mat2 <- t(pctmat.cutoff.E2toE1[-nrow(pctmat.cutoff.E2toE1),]) #so transpose E2toE1
  mat.bi <- mat1+mat2
  ## unassigned row
  mat.bi <- rbind(mat.bi, 2*as.numeric(colSums(mat.bi)==0))
  ## rename colnames and rownames
  rownames(mat.bi) <- gsub("ref",name.E2,rownames(pctmat.cutoff.E1toE2))
  colnames(mat.bi) <- gsub("query",name.E1,colnames(pctmat.cutoff.E1toE2))

  ## plot
  if(is.null(main)) main <- "scmap"
  if(reorder) mat.bi <- reorder(mat.bi)
  pheatmap(mat.bi,
           color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(1,1,3,7)]))(4),
           breaks=seq(0,2,length.out=4),
           legend_breaks=0:2, legend_labels=c("No match", "One-way match", "Two-way Match"),
           cluster_rows=F, cluster_cols=F,
           gaps_row=nrow(mat.bi)-1,
           cellwidth=cellwidth, cellheight=cellheight,
           main=main, filename=filename,
           ...)
  cat("pct.cutoff = ", pct.cutoff, "\n")
  if(return.value) return(mat.bi)
}


######################################################################################################################
### Seurat ###########################################################################################################
######################################################################################################################


#######################
## match_by_seurat() ##
#######################

match_by_seurat <- function(sce.query, sce.ref, #imputation=FALSE,
                            filter.size=10, filter.fscore=NULL, #filter clusters
                            use.markers=FALSE, dims=1:30){

  ## filtering small or low fscore clusters
  cat("Filtering small clusters: clusters with less than", filter.size, "cells are not considered. \n")
  if(!is.null(filter.fscore)) cat("Filtering low F-beta score clusters in the reference experiment only. \n")
  sce.ref <- filter_cluster(sce.ref, filter.size=filter.size, filter.fscore=filter.fscore)
  sce.query <- filter_cluster(sce.query, filter.size=filter.size, filter.fscore=NULL)

  ## imputation
  # if(imputation){
  #   sce.ref <- impute_dropout(sce.ref)
  #   cat("Imputation done. \n")
  # }

  ## extract info from sce.objects
  oo.query <- sce.query@metadata$cluster_order
  oo.ref <- sce.ref@metadata$cluster_order

  ## make data into seurat objects in a list
  cat("Preparing data for Seurat... \n")
  mat <- merge(assay(sce.query), assay(sce.ref), by="row.names")
  rownames(mat) <- mat$Row.names
  mat <- as.matrix(mat[,-1])
  data <- as(mat, "dgCMatrix")

  metadata <- data.frame("subsample"=rep(c("query","ref"), c(ncol(sce.query),ncol(sce.ref))),
                         "celltype"=c(as.character(colData(sce.query)$cluster_membership),
                                      as.character(colData(sce.ref)$cluster_membership)),
                         stringsAsFactors=FALSE)
  rownames(metadata) <- colnames(mat)

  seurat.object <- Seurat::CreateSeuratObject(counts=data, meta.data=metadata)
  seurat.object.list <- Seurat::SplitObject(object=seurat.object, split.by="subsample")

  ## manually select features of NS-Forest marker genes
  if(use.markers){
    seurat.object.list[["query"]]@assays$RNA@var.features <- unique(sce.query@metadata$cluster_marker_info$markerGene)
    seurat.object.list[["ref"]]@assays$RNA@var.features <- unique(sce.ref@metadata$cluster_marker_info$markerGene)
    cat("Customized markers are used. \n")
  }
  ## select variable features
  else for (i in 1:length(x=seurat.object.list)) {
    seurat.object.list[[i]] <- Seurat::FindVariableFeatures(object=seurat.object.list[[i]],
                                                            selection.method="vst", nfeatures=2000, verbose=FALSE)
  }

  ### matching ###
  cat("Matching by Seurat... \n")
  object.ref <- seurat.object.list[[c("ref")]]
  object.query <- seurat.object.list[[c("query")]]
  ## PCA used by default in the following step, which is recommended for scRNA-seq
  anchors <- Seurat::FindTransferAnchors(reference=object.ref, query=object.query, dims=dims)
  predictions <- Seurat::TransferData(anchorset=anchors, refdata=object.ref$celltype, dims=dims)

  ## calculate Seurat scores averaged by cluster
  df.predictions <- predictions %>% rownames_to_column()
  df.predictions$rowname <- gsub(".x","",df.predictions$rowname)
  df.membership <- colData(sce.query) %>% as.data.frame() %>% rownames_to_column()
  df.scores <- df.predictions %>% left_join(df.membership) %>% group_by(cluster_membership) %>%
    select(-c(rowname, predicted.id, prediction.score.max)) %>%
    summarise_all(mean) # average scores by cluster
  mymat.scores <- df.scores %>% select(starts_with("prediction")) %>% as.matrix() %>% t()
  colnames(mymat.scores) <- paste0("query.", df.scores$cluster_membership)
  scores <- mymat.scores[paste0("prediction.score.",oo.ref),paste0("query.",oo.query)]

  ## prediction.celltype
  object.query <- Seurat::AddMetaData(object=object.query, metadata=predictions)
  prediction.celltype <- cbind(object.query$predicted.id, object.query$celltype)
  colnames(prediction.celltype) <- c("prediction", "celltype")

  ## match.lst
  celltypes <- unique(prediction.celltype[,"celltype"])
  match.lst <- vector("list",length=length(celltypes))
  names(match.lst) <- celltypes
  for(celltype in celltypes){
    ind <- prediction.celltype[,"celltype"]==celltype
    match.lst[[celltype]] <- table(prediction.celltype[ind,1])
  }
  match.lst <- match.lst[oo.query]

  cluster.sizes.query <- table(colData(sce.query)$cluster_membership)[names(match.lst)]
  pctmat <- match2mat.seurat(match.lst, oo.ref, cluster.sizes.query)

  ## all results
  all.results <- data.frame("query.cell"=rownames(prediction.celltype), "query.cluster"=prediction.celltype[,"celltype"],
                            "predicted.ref.cluster"=prediction.celltype[,"prediction"])

  ## return
  return(list("method"="Seurat v3 Integration - Standard Workflow",
              "pctmat"=pctmat, "prediction.scores"=scores, "all.results"=all.results))
}

##--------------------##
## match2mat.seurat() ##
##--------------------##

match2mat.seurat <- function(match.lst, oo.ref, cluster.sizes.query){
  df.match <- rownames_to_column(data.frame("freq"=unlist(match.lst)))
  df.all <- merge(names(match.lst), oo.ref, all=TRUE) %>%
    mutate(query.ref=paste0(x,".",y)) %>%
    left_join(df.match, by=c("query.ref"="rowname")) %>%
    mutate(freq=replace_na(freq,0))
  mat <- matrix(df.all$freq, nrow=length(oo.ref), byrow=TRUE)

  unassigned <- cluster.sizes.query - colSums(mat)
  mat <- rbind(mat,unassigned)

  rownames(mat) <- c(paste0("ref.",oo.ref), "unassigned") #row=ref
  colnames(mat) <- paste0("query.",names(match.lst)) #col=query

  mat.pct <- sweep(mat,2,cluster.sizes.query,"/")
  return(mat.pct)
}

###################
## plot_seurat() ##
###################

plot_seurat <- function(rst.seurat, type="matches", pct.cutoff=0.3,
                        reorder=TRUE, return.value=FALSE,
                        cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){
  pctmat <- rst.seurat$pctmat
  pctmat.cutoff <- cutoff.pctmat(pctmat, pct.cutoff=pct.cutoff)
  if(reorder){
    pctmat.cutoff <- reorder(pctmat.cutoff)
    pctmat <- pctmat[,colnames(pctmat.cutoff)]
  }
  if(type=="matches"){
    if(is.null(main)) main <- "Seurat"
    ## heatmap
    pheatmap(pctmat.cutoff,
             color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(3,3,7)]))(3),
             breaks=seq(0,1,length.out=3),
             legend_breaks=c(0,1), legend_labels=c("No match", "Match"),
             cluster_rows=F, cluster_cols=F,
             gaps_row=nrow(pctmat.cutoff)-1,
             cellwidth=cellwidth, cellheight=cellheight,
             main=main, filename=filename,
             ...)
    cat("pct.cutoff = ", pct.cutoff, "\n")
    ## output
    if(return.value) return(pctmat.cutoff)
  }
  if(type=="pctmat"){
    # if(is.null(main)) main <- "Distribution of % matched cells"
    pctmat <- pctmat[-nrow(pctmat),]
    df <- tibble(pct=as.vector(pctmat),
                 query_cluster = rep(colnames(pctmat), each=nrow(pctmat))) %>%
      mutate(query_cluster = fct_relevel(query_cluster, colnames(pctmat)))
    g <- ggplot(df, aes(x=query_cluster, y=pct)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 270, hjust = 0)) +
      geom_hline(linetype = "dashed", yintercept = pct.cutoff, color = "red") +
      scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      xlab("Query cluster") + ylab("% matched cells")
    ## save plot or plot on device
    if(!is.na(filename)) ggsave(filename, g, width=ncol(pctmat)*.2, height=5)
    else plot(g)
    cat("pct.cutoff = ", pct.cutoff, "\n")
    ## output
    if(return.value) return(pctmat)
  }
  if(type=="prediction.scores"){
    prediction.scores <- rst.seurat$prediction.scores
    if(reorder) prediction.scores <- prediction.scores[,colnames(pctmat.cutoff)]
    if(is.null(main)) main <- "Seurat prediction score"
    pheatmap(prediction.scores,
             color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(101),
             breaks=seq(0,1,length.out=101),
             cluster_rows=F, cluster_cols=F,
             cellwidth=cellwidth, cellheight=cellheight,
             main=main, filename=filename,
             ...)
    if(return.value) return(prediction.scores)
  }
}

######################
## plot_bi_seurat() ##
######################

plot_bi_seurat <- function(rst.seurat.E1toE2, rst.seurat.E2toE1, name.E1="E1.", name.E2="E2.",
                           pct.cutoff=0.3,
                           reorder=TRUE, return.value=FALSE,
                           cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){
  ## get binary matrices for plotting
  pctmat.cutoff.E1toE2 <- cutoff.pctmat(rst.seurat.E1toE2$pctmat, pct.cutoff=pct.cutoff)
  pctmat.cutoff.E2toE1 <- cutoff.pctmat(rst.seurat.E2toE1$pctmat, pct.cutoff=pct.cutoff)

  ## combine two matrices to one two-way matrix
  mat1 <- pctmat.cutoff.E1toE2[-nrow(pctmat.cutoff.E1toE2),] #use E1toE2 as the framework for final plot
  mat2 <- t(pctmat.cutoff.E2toE1[-nrow(pctmat.cutoff.E2toE1),]) #so transpose E2toE1
  mat.bi <- mat1+mat2
  ## unassigned row
  mat.bi <- rbind(mat.bi, 2*as.numeric(colSums(mat.bi)==0))
  ## rename colnames and rownames
  rownames(mat.bi) <- gsub("ref.", name.E2, rownames(pctmat.cutoff.E1toE2))
  colnames(mat.bi) <- gsub("query.", name.E1, colnames(pctmat.cutoff.E1toE2))

  ## plot
  if(is.null(main)) main <- "Seurat"
  if(reorder) mat.bi <- reorder(mat.bi)
  pheatmap(mat.bi,
           color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(1,1,3,7)]))(4),
           breaks=seq(0,2,length.out=4),
           legend_breaks=0:2, legend_labels=c("No match", "One-way match", "Two-way Match"),
           cluster_rows=F, cluster_cols=F,
           gaps_row=nrow(mat.bi)-1,
           cellwidth=cellwidth, cellheight=cellheight,
           main=main, filename=filename,
           ...)
  cat("pct.cutoff = ", pct.cutoff, "\n")
  if(return.value) return(mat.bi)
}

plot_seurat_cell2cluster <- function(rst, type="match.prop",
                                     reorder=TRUE, return.value=FALSE,
                                     main=NULL, filename=NA, width=NULL, height=NULL){
  if(type=="match.prop"){
    tab.match <- table(rst$all.results$predicted.ref.cluster, rst$all.results$query.cluster) #query subclass in columns, prediction in rows
    tab.match.prop <- sweep(tab.match,2,colSums(tab.match),"/") #column sums should be 1

    ## order reference cluster orders in rows
    clusterNames.ref <- gsub("ref.","",rownames(rst$pctmat)) #reference cluster names IN ORDER
    oo <- match(clusterNames.ref, rownames(tab.match.prop))
    oo.names <- rownames(tab.match.prop)[oo] %>% na.omit() #some ref clusters may not have matched query cells
    tab.match.prop <- tab.match.prop[oo.names,]

    ## reorder columns for query
    if(reorder) tab.match.prop <- reorder(tab.match.prop)

    ## for ggplot
    long.tab.match.prop <- tab.match.prop %>% as.data.frame() %>%
      select(query.cluster=Var2, match=Var1, Prop=Freq) %>%
      mutate(match=factor(match, levels = rev(clusterNames.ref)))

    ## plot
    if(is.null(main)) main <- "Seurat cell-to-cluster"
    g <- ggplot(long.tab.match.prop, aes(x=query.cluster, y=match, size=Prop, fill=Prop)) +
      geom_point(alpha=0.7, shape=21, color="black") +
      scale_size_continuous(range = c(0, 10), limits = c(0, 1)) +
      scale_fill_viridis(option="D", guide = "legend", limits = c(0, 1)) +
      scale_y_discrete(drop=FALSE) + #show all ref clusters even if no match
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      ggtitle(main)

    ## save plot or plot on device
    if(is.null(width)) width <- ncol(tab.match.prop)*.2+.5
    if(is.null(height)) height <- nrow(tab.match.prop)*.2
    if(!is.na(filename)) ggsave(filename, g, width=width, height=height)
    else plot(g)

    ## output
    if(return.value){
      return(list("plotted.values"=tab.match.prop))
    }
  }
}

