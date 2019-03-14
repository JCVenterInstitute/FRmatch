

#' Plotting function for FRmatch results
#'
#' This function facilitates visualization of \code{FRmatch} results.
#'
#' @param rst.FRmatch The output retults from the \code{\link[FRmath]{FRmatch}} function.
#' @param type Character variable that specifies which element from the output list to be ploted.
#' Default: \code{"matches"}, which returns the recommended matches found by FRmatch p-values thresholded at \code{sig.level}.
#' Other options: \code{"pmat", "statmat"}.
#' @param sig.level Numeric variable that sets the significance level for FRmatch p-values, above which is a match.
#' Default value: \code{0.05}.
#' @param cellwidth,cellheight,... Additional plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @seealso \link[pheatmap]{pheatmap}.
#' @export

plot_FRmatch <- function(rst.FRmatch, type="matches", sig.level=0.05, cellwidth=10, cellheight=10, ...){
  if(type=="matches"){
    pmat.cutoff <- myfun.cutoff(rst.FRmatch$pmat, cutoff=sig.level)
    pheatmap::pheatmap(pmat.cutoff,
                       color = colorRampPalette(rev(RColorBrewer::brewer.pal(n=7, name="RdYlBu")[c(1,1,7)]))(3),
                       breaks = seq(0,1,length.out=3),
                       legend_breaks=c(0,1), legend_labels=c("Not match", "Match"),
                       cluster_rows=F, cluster_cols=F,
                       gaps_row=nrow(pmat.cutoff)-1,
                       cellwidth=cellwidth, cellheight=cellheight,
                       main="Recommended matches",
                       ...)
  }
  if(type=="pmat"){
    pheatmap::pheatmap(rst.FRmatch$pmat,
             cluster_rows=F, cluster_cols=F,
             cellwidth=cellwidth, cellheight=cellheight,
             main="P-values",
             ...)
  }
  if(type=="statmat"){
    pheatmap::pheatmap(rst.FRmatch$statmat,
             cluster_rows=F, cluster_cols=F,
             cellwidth=cellwidth, cellheight=cellheight,
             main="FR statistics",
             ...)
  }
}


myfun.cutoff <- function(pmat, cutoff){
  out <- matrix(as.numeric(pmat>cutoff),nrow=nrow(pmat))
  out <- rbind(out, as.numeric(colSums(out)==0))
  colnames(out) <- colnames(pmat)
  rownames(out) <- c(rownames(pmat), "unassigned")
  return(out)
}
