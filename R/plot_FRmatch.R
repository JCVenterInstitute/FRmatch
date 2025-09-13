
#' Plotting function for FR-Match results
#'
#' This function takes in the \code{\link[FRmatch]{FRmatch}} output and generates plots for the matching results.
#'
#' @param rst.FRmatch The \code{\link[FRmatch]{FRmatch}} output.
#' @param type If \code{type="matches"} (default), it plots the one-way matches.
#' If \code{type="padj"}, it plots the distribution of adjusted p-values.
#' @param p.adj.method P-value adjustment method for multiple hypothesis testing correction. Default: \code{"BY"}.
#' For more options, please see \code{\link[stats]{p.adjust.methods}}.
#' @param sig.level Numeric variable that specifies the significance level of adjusted p-value. A MATCH is >\code{sig.level}.
#' Default: \code{0.05}.
#' @param reorder Boolean variable indicating if to reorder the columns so that matches are aligned along the diagonal.
#' It improves the interpretability of the one-way match plot. Default: \code{TRUE}.
#' @param ignore.unassigned Boolean variable indicating if to skip the columns of unassigned query clusters
#' in the \code{type="matches"} plot. Default: \code{FALSE}. If \code{TRUE}, number of ignored columns reported in the unassigned row.
#' @param return.value Boolean variable indicating if to return the plotted values. Default: \code{FALSE}.
#' @param cellwidth,cellheight,main,filename,... Plotting parameters passed to \code{\link[graphics]{hist}}.
#'
#' @return If \code{return.value = TRUE}, a matrix of one-way matching values 1 = match, and 0 = no match, or a matrix of adjusted p-values.
#'
#' @seealso \code{\link[FRmatch]{plot_bi_FRmatch}}.
#'
#' @importFrom forcats fct_relevel
#'
#' @export

plot_FRmatch <- function(rst.FRmatch, type="matches", p.adj.method="BY", sig.level=0.05,
                         reorder=TRUE, ignore.unassigned=FALSE, return.value=FALSE,
                         cellwidth=10, cellheight=10, main=NULL, filename=NA, ...){
  ## calculate adjusted p-values and determine matches
  pmat.adj <- padj.FRmatch(rst.FRmatch$pmat, p.adj.method=p.adj.method)
  pmat.cutoff <- cutoff.FRmatch(rst.FRmatch$pmat, p.adj.method=p.adj.method, sig.level=sig.level)

  ## reorder
  if(reorder){
    pmat.cutoff <- reorder(pmat.cutoff)
    pmat.adj <- pmat.adj[,colnames(pmat.cutoff)]
  }

  ## plot
  if(type=="matches"){
    if(is.null(main)) main <- "FR-Match cluster-to-cluster"
    if(ignore.unassigned){
      ind <- pmat.cutoff["unassigned",]==1 #unassigned columns
      pmat.cutoff <- pmat.cutoff[,!ind]
      rownames(pmat.cutoff)[nrow(pmat.cutoff)] <- paste0("unassigned (", sum(ind),")")
    }
    ## heatmap
    pheatmap(pmat.cutoff,
             color = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")[c(3,3,7)]))(3),
             breaks = seq(0,1,length.out=3),
             legend_breaks=c(0,1), legend_labels=c("No match", "Match"),
             cluster_rows=F, cluster_cols=F,
             gaps_row=nrow(pmat.cutoff)-1,
             cellwidth=cellwidth, cellheight=cellheight,
             main=main,
             filename=filename,
             ...)
    ## output
    if(return.value) return(pmat.cutoff)
  }

  ## plot
  if(type=="padj"){
    df <- tibble(padj=as.vector(pmat.adj),
                 query.cluster = rep(colnames(pmat.adj), each=nrow(pmat.adj))) %>%
      mutate(query.cluster = fct_relevel(query.cluster, colnames(pmat.adj)))
    g <- ggplot(df, aes(x=query.cluster, y=padj)) +
      geom_boxplot() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 270, hjust = 0)) +
      geom_hline(linetype = "dashed", yintercept = sig.level, color = "red") +
      scale_y_continuous(breaks=seq(0,1,0.2), limits=c(0,1)) +
      xlab("Query cluster") + ylab("Adjusted p-value")
    ## save plot or plot on device
    if(!is.na(filename)) ggsave(filename, g, width=ncol(pmat.adj)*.2, height=5)
    else plot(g)
    ## output
    if(return.value) return(pmat.adj)
  }
}


##------below are some utility functions for the main function---------------------------------


######################
## cutoff.FRmatch() ##
######################
## This function determines matches by cutting off the p-values from FRmatch, and
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




