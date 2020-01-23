
#' Bilateral matching plot of FR-Match results
#'
#' This function combines two sets of reciprocal \code{\link[FRmath]{FRmatch}} outputs (i.e. experiment 1 (hereinafter, E1) as the query
#' dataset mapping to experiment 2 (hereinafter, E2) as the reference dataset, and E2 mapping to E1), and plots the two-way matches,
#' one-way matches, and no matches determined by FR-Match.
#'
#' @param rst.FRmatch.E1toE2,rst.FRmatch.E2toE1 The \code{\link[FRmath]{FRmatch}} outputs.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param p.adj.method See \code{\link[FRmath]{plot_FRmatch}}.
#' @param sig.level See \code{\link[FRmath]{plot_FRmatch}}.
#' @param reorder See \code{\link[FRmath]{plot_FRmatch}}.
#' @param return.value Logical variable indicating if to return plotted values. Default: \code{FALSE}.
#' @param cellwidth,cellheight,main,filename,... Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return If \code{return.value = TRUE}, a matrix of matching results. 2 = two-way match, 1 = one-way match, and 0 = no match.
#'
#' @seealso \code{\link[FRmatch]{plot_FRmatch}}.
#'
#' @import pheatmap
#' @importFrom RColorBrewer brewer.pal
#'
#' @export

plot_bilateral_FRmatch <- function(rst.FRmatch.E1toE2, rst.FRmatch.E2toE1, name.E1="E1", name.E2="E2",
                                   p.adj.method="BY", sig.level=0.05,
                                   reorder=TRUE, return.value=FALSE,
                                   cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){

  ## get binary matrices for plotting
  pmat.cutoff.E1toE2 <- cutoff.FRmatch(rst.FRmatch.E1toE2$pmat, p.adj.method=p.adj.method, sig.level=sig.level)
  pmat.cutoff.E2toE1 <- cutoff.FRmatch(rst.FRmatch.E2toE1$pmat, p.adj.method=p.adj.method, sig.level=sig.level)

  ## combine two matrices to one bilateral matrix
  mat1 <- pmat.cutoff.E1toE2[-nrow(pmat.cutoff.E1toE2),] #use E1toE2 as the framework for final plot
  mat2 <- t(pmat.cutoff.E2toE1[-nrow(pmat.cutoff.E2toE1),]) #so transpose E2toE1
  mat.bi <- mat1+mat2
  ## unassigned row
  mat.bi <- rbind(mat.bi, 2*as.numeric(colSums(mat.bi)==0))
  ## rename colnames and rownames
  rownames(mat.bi) <- gsub("ref",name.E2,rownames(pmat.cutoff.E1toE2))
  colnames(mat.bi) <- gsub("query",name.E1,colnames(pmat.cutoff.E1toE2))

  ## plot
  if(is.null(main)) main <- "Bilatreral matches"
  if(reorder) mat.bi <- reorder(mat.bi)
  pheatmap(mat.bi,
           color = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(1,1,3,7)]))(4),
           breaks = seq(0,2,length.out=4),
           legend_breaks=0:2, legend_labels=c("No match", "One-way match", "Two-way Match"),
           cluster_rows=F, cluster_cols=F,
           gaps_row=nrow(mat.bi)-1,
           cellwidth=cellwidth, cellheight=cellheight,
           main=main,
           filename=filename,
           ...)

  ## output
  if(return.value) return(mat.bi)
}


