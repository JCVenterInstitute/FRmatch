
#' Visualization of bilateral matching restuls of FRmatch
#'
#' This function combines two reciprocal \code{FRmatch} result objects (i.e. experiment 1 (hereinafter, E1) as query
#' mapping to E2 as reference, and E2 as query mapping to E1 as reference), and plots the matches that are identified
#' in both ways or just in one way or no match.
#'
#' @param rst.FRmatch.E1toE2,rst.FRmatch.E2toE1 Retults outputed from the \code{\link[FRmath]{FRmatch}} function.
#' @param prefix.E1,prefix.E2 Prefix names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param p.adj.method P-value adjustment method for multiple comparison correction. Default: \code{"BY"}.
#' Please see \code{\link[stats]{p.adjust.methods}}.
#' @param sig.level Numeric variable that sets the significance level for the adjusted p-values, above which is a match.
#' Default value: \code{0.05}.
#' @param reorder Logical variable indicating if to reorder the columns of heatmap for better interpretability
#' (aligning matches in the diagonal). Default: \code{TRUE}.
#' @param return.value Logical variable indicating if to return the values correponding to the plot. Default: \code{FALSE}.
#' @param cellwidth,cellheight,main,... Additional plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return Optionally, a numeric matrix corresponding to the values on the plot.
#'
#' @seealso \link[FRmatch]{plot_FRmatch}.
#' @export

plot_bilateral_FRmatch <- function(rst.FRmatch.E1toE2, rst.FRmatch.E2toE1, prefix.E1="E1", prefix.E2="E2",
                                   p.adj.method="BY", sig.level=0.05,
                                   reorder=TRUE, return.value=FALSE,
                                   cellwidth=10, cellheight=10, main=NULL, ...){
  ## get binary matrices for plotting
  pmat.cutoff.E1toE2 <- FRmatch:::cutoff.FRmatch(rst.FRmatch.E1toE2$pmat, p.adj.method=p.adj.method, sig.level=sig.level)
  pmat.cutoff.E2toE1 <- FRmatch:::cutoff.FRmatch(rst.FRmatch.E2toE1$pmat, p.adj.method=p.adj.method, sig.level=sig.level)

  ## combine two matrices to one bilateral matrix
  mat1 <- pmat.cutoff.E1toE2[-nrow(pmat.cutoff.E1toE2),] #use E1toE2 as the framework for final plot
  mat2 <- t(pmat.cutoff.E2toE1[-nrow(pmat.cutoff.E2toE1),]) #so transpose E2toE1
  mat.bi <- mat1+mat2
  ## unassigned row
  mat.bi <- rbind(mat.bi, 2*as.numeric(colSums(mat.bi)==0))
  ## rename colnames and rownames
  rownames(mat.bi) <- gsub("ref",prefix.E2,rownames(pmat.cutoff.E1toE2))
  colnames(mat.bi) <- gsub("query",prefix.E1,colnames(pmat.cutoff.E1toE2))

  ## plot
  if(is.null(main)) main <- "Bilatreral matches"
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


