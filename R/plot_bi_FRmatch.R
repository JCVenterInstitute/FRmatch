
#' Plotting function for bi-directional FR-Match results
#'
#' This function combines and plots two sets of reciprocal \code{\link[FRmatch]{FRmatch}} results (i.e. experiment 1 (E1) query
#' to experiment 2 (E2) reference with E2 markers, and E2 query to E1 reference with E1 markers).
#'
#' @param rst.FRmatch.E1toE2,rst.FRmatch.E2toE1 The \code{\link[FRmatch]{FRmatch}} outputs.
#' @param prefix Prefix names for query and reference clusters used in \code{\link[FRmatch]{FRmatch}}.
#' @param name.E1,name.E2 Customized names with delimiter for E1 and E2 to be used in this figure. Default: \code{"E1."} and \code{"E2."}, respectively.
#' @param p.adj.method See \code{\link[FRmatch]{plot_FRmatch}}.
#' @param sig.level See \code{\link[FRmatch]{plot_FRmatch}}.
#' @param reorder See \code{\link[FRmatch]{plot_FRmatch}}.
#' @param two.way.only Boolean variable indicating if to plot two-way matches only. Default: \code{FALSE}.
#' @param return.value Boolean variable indicating if to return the plotted values. Default: \code{FALSE}.
#' @param cellwidth,cellheight,main,filename,... Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return If \code{return.value = TRUE}, a matrix of two-way matching values 2 = two-way match, 1 = one-way match, and 0 = no match.
#'
#' @seealso \code{\link[FRmatch]{plot_FRmatch}}.
#'
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#'
#' @export

plot_bi_FRmatch <- function(rst.FRmatch.E1toE2, rst.FRmatch.E2toE1,
                            prefix=c("query.","ref."), name.E1="E1.", name.E2="E2.",
                            p.adj.method="BY", sig.level=0.05,
                            reorder=TRUE, two.way.only=FALSE, return.value=FALSE,
                            cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){

  ## get binary matrices for plotting
  pmat.cutoff.E1toE2 <- cutoff.FRmatch(rst.FRmatch.E1toE2$pmat, p.adj.method=p.adj.method, sig.level=sig.level)
  pmat.cutoff.E2toE1 <- cutoff.FRmatch(rst.FRmatch.E2toE1$pmat, p.adj.method=p.adj.method, sig.level=sig.level)

  ## combine two matrices to one two-way matrix
  mat1 <- pmat.cutoff.E1toE2[-nrow(pmat.cutoff.E1toE2),] #use E1toE2 as the framework for final plot
  mat2 <- t(pmat.cutoff.E2toE1[-nrow(pmat.cutoff.E2toE1),]) #so transpose E2toE1
  mat.bi <- mat1+mat2
  if(two.way.only) mat.bi <- matrix(2*as.numeric(mat.bi==2), nrow(mat1), ncol(mat1))
  ## unassigned row
  mat.bi <- rbind(mat.bi, 2*as.numeric(colSums(mat.bi)==0))
  ## rename colnames and rownames
  rownames(mat.bi) <- gsub(prefix[2],name.E2,rownames(pmat.cutoff.E1toE2))
  colnames(mat.bi) <- gsub(prefix[1],name.E1,colnames(pmat.cutoff.E1toE2))

  ## plot
  if(is.null(main)) main <- "FR-Match cluster-to-cluster"
  if(reorder) mat.bi <- reorder(mat.bi)
  if(two.way.only){
    pheatmap(mat.bi,
             color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(1,1,7)]))(3),
             breaks=seq(0,2,length.out=3),
             legend_breaks=c(0,2),
             legend_labels=c("No match", "Two-way match"),
             cluster_rows=F, cluster_cols=F,
             gaps_row=nrow(mat.bi)-1,
             cellwidth=cellwidth, cellheight=cellheight,
             main=main,
             filename=filename,
             ...)
  } else
    pheatmap(mat.bi,
             color=colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(1,1,3,7)]))(4),
             breaks=seq(0,2,length.out=4),
             legend_breaks=0:2, legend_labels=c("No match", "One-way match", "Two-way match"),
             cluster_rows=F, cluster_cols=F,
             gaps_row=nrow(mat.bi)-1,
             cellwidth=cellwidth, cellheight=cellheight,
             main=main,
             filename=filename,
             ...)

  ## output
  if(return.value) return(mat.bi)
}


