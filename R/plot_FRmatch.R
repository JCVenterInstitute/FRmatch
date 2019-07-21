
#' Visulization of FRmatch result
#'
#' This function plots results outputed from the \code{\link[FRmath]{FRmatch}} function. By default (\code{type="matches"}),
#' it returns the recomended matches determined by adjusted p-values using \code{p.adj.method} on FRmatch p-values at the
#' significant level of \code{sig.level}.
#'
#' @param rst.FRmatch The result object outputed from the \code{\link[FRmath]{FRmatch}} function.
#' @param type Character variable that specifies which element from the output list to be ploted.
#' Default: \code{"matches"}, which returns the recommended matches identified by FRmatch using \code{p.adj.method}
#' and \code{sig.level}. Other options: \code{"pmat", "statmat"}.
#' @param p.adj.method P-value adjustment method for multiple comparison correction. Default: \code{"BY"}.
#' Please see \code{\link[stats]{p.adjust.methods}}.
#' @param sig.level Numeric variable that sets the significance level for the adjusted p-values, above which is a match.
#' Default value: \code{0.05}.
#' @param reorder Logical variable indicating if to reorder the columns of heatmap for better interpretability
#' (aligning matches in the diagonal). This parameter only applies to \code{type="matches"}, default: \code{TRUE}.
#' @param return.value Logical variable indicating if to return the values correponding to the "matches" plot.
#' This parameter only applies to \code{type="matches"}, default: \code{FALSE}.
#' @param cellwidth,cellheight,main,... Additional plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return Optionally, a numeric matrix corresponding to the values on the plot.
#'
#' @seealso \link[FRmatch]{FRmatch}.
#' @importFrom dplyr %>%
#' @export

####################
## plot_FRmatch() ##
####################

plot_FRmatch <- function(rst.FRmatch, type="matches", p.adj.method="BY", sig.level=0.05,
                         reorder=TRUE, return.value=FALSE,
                         cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){
  pmat.adj <- padj.FRmatch(rst.FRmatch$pmat, p.adj.method=p.adj.method)
  pmat.cutoff <- cutoff.FRmatch(rst.FRmatch$pmat, p.adj.method=p.adj.method, sig.level=sig.level)
  if(reorder){
    pmat.cutoff <- reorder(pmat.cutoff)
    pmat.adj <- pmat.adj[,colnames(pmat.cutoff)]
  }
  if(type=="matches"){
    if(is.null(main)) main <- "Recommended matches"
    pheatmap::pheatmap(pmat.cutoff,
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")[c(3,3,7)]))(3),
                       breaks = seq(0,1,length.out=3),
                       legend_breaks=c(0,1), legend_labels=c("No match", "Match"),
                       cluster_rows=F, cluster_cols=F,
                       gaps_row=nrow(pmat.cutoff)-1,
                       cellwidth=cellwidth, cellheight=cellheight,
                       main=main,
                       filename=filename,
                       ...)
  }
  if(type=="padj"){
    # if(is.null(main)) main <- "Distribution of adjusted p-values"
    df <- tibble::tibble(padj=as.vector(pmat.adj),
                         query_cluster = rep(colnames(pmat.adj), each=nrow(pmat.adj))) %>%
      dplyr::mutate(query_cluster = forcats::fct_relevel(query_cluster, colnames(pmat.adj)))
    g <- ggplot2::ggplot(df, ggplot2::aes(x=query_cluster, y=padj)) +
      ggplot2::geom_boxplot() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
      ggplot2::geom_hline(linetype = "dashed", yintercept = sig.level, color = "red")
    plot(g)
    if(!is.na(filename)) ggplot2::ggsave(filename, g, width=ncol(pmat.adj)*.2, height=5)
  }
  if(return.value){
    return(list("matches"=pmat.cutoff, "padj"=pmat.adj))
  }
  # if(type=="pmat"){
  #   if(is.null(main)) main <- "P-values"
  #   pheatmap::pheatmap(rst.FRmatch$pmat,
  #            cluster_rows=F, cluster_cols=F,
  #            cellwidth=cellwidth, cellheight=cellheight,
  #            main=main,
  #            ...)
  # }
  # if(type=="pmat.rk"){
  #   if(is.null(main)) main <- "Rank of p-values"
  #   pheatmap::pheatmap(apply(rst.FRmatch$pmat, 2, rank),
  #                      cluster_rows=F, cluster_cols=F,
  #                      cellwidth=cellwidth, cellheight=cellheight,
  #                      main=main,
  #                      ...)
  # }
  # if(type=="statmat"){
  #   if(is.null(main)) main <- "FR statistics"
  #   pheatmap::pheatmap(rst.FRmatch$statmat,
  #            cluster_rows=F, cluster_cols=F,
  #            cellwidth=cellwidth, cellheight=cellheight,
  #            main=main,
  #            ...)
  # }
}

##------below are some utility function for the main function---------------------------------


######################
## cutoff.FRmatch() ##
######################
## This function determines matches by cuting off the p-values from FRmatch, and
## adds the "unassigned" row in the bottom

cutoff.FRmatch <- function(pmat, p.adj.method, sig.level){
  pmat.adj <- padj.FRmatch(pmat=pmat, p.adj.method=p.adj.method)
  out <- matrix(as.numeric(pmat.adj>=sig.level), nrow=nrow(pmat.adj))
  out <- rbind(out, as.numeric(colSums(out)==0))
  colnames(out) <- colnames(pmat)
  rownames(out) <- c(rownames(pmat), "unassigned")
  return(out)
}

padj.FRmatch <- function(pmat, p.adj.method){
  pvals.adj <- p.adjust(pmat, method=p.adj.method)
  pmat.adj <- matrix(pvals.adj, nrow=nrow(pmat))
  colnames(pmat.adj) <- colnames(pmat)
  rownames(pmat.adj) <- c(rownames(pmat))
  return(pmat.adj)
}

###############
## reorder() ##
###############
## This function reorders the heatmap for better interpretability (diagonalizing the matches)
## Column order: if we could start with the one-to-one and then many-to-one and then unassigned

reorder <- function(MATRIX){
  tmp <- colSums((MATRIX>0)[-nrow(MATRIX),]) #last row is "unassigned" - don't count't
  tmp[tmp>0] <- 1
  tmp[tmp==0] <- max(tmp)+1
  ord <- order(tmp, apply(MATRIX[-nrow(MATRIX),],2,which.max))
  return(MATRIX[,ord])
}




