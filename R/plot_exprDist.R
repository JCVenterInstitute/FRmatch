
#' Gene expression data distribution plot
#'
#' This function plots the expression data distributions of the two single cell datasets (e.g. query and reference) to be compared.
#'
#' @param sce.E1,sce.E2 Data objects, namely E1 and E2.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param breaks,xlim,ylim Plotting parameters passed to histogram plot.
#' @param filename File name if to save the plot. Default: \code{NA}, not to save the plot.
#' @param width,height Width and height for saved plot.
#' @param ... Additional arguments passed to \code{\link[graphics]{hist}}
#'
#' @export

plot_exprDist <- function(sce.E1, sce.E2, markers=NULL, name.E1="E1", name.E2="E2",
                          breaks=20, xlim=c(0,10), ylim=c(0,2),
                          filename=NA, width=10, height=5, ...){
  ## to save pdf
  if(!is.na(filename)){pdf(filename, width=width, height=height)}

  ## get common markers for plot
  if(is.null(markers)){
    markers.common <- intersect(rownames(sce.E1), rownames(sce.E2))
  } else
    markers.common <- intersect(intersect(markers, rownames(sce.E1)), rownames(sce.E2))

  ## plot
  par(mfrow=c(1,2), mar=c(3,4,3,2))
  hist(logcounts(sce.E1)[markers.common,], freq=F, xlab="",
       breaks=breaks, xlim=xlim, ylim=ylim, main=name.E1, ...)
  ss <- summary(as.vector(logcounts(sce.E1)))
  legend("topright", paste(names(ss),"=", round(ss,3)), bty="n")
  hist(logcounts(sce.E2)[markers.common,], freq=F, xlab="",
       breaks=breaks, xlim=xlim, ylim=ylim, main=name.E2, ...)
  ss <- summary(as.vector(logcounts(sce.E2)))
  legend("topright", paste(names(ss),"=", round(ss,3)), bty="n")

  ## to close pdf
  if(!is.na(filename)){dev.off()}
}
