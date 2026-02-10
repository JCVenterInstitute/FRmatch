
#' Plotting function for FR-Match cell2cluster results
#'
#' This function takes in the \code{\link[FRmatch]{FRmatch_cell2cluster}} output and generates plots for the matching results.
#'
#' @param rst.cell2cluster The \code{\link[FRmatch]{FRmatch_cell2cluster}} output.
#' @param type If \code{type="match"} (default), it plots the proportion of cells in the query cluster matched to the reference cluster.
#' If \code{type=="padj"}, it plots the distribution of adjusted p-values.
#' @param p.adj.method P-value adjustment method for multiple hypothesis testing correction. Default: \code{"BH"}.
#' For more options, please see \code{\link[stats]{p.adjust.methods}}.
#' @param sig.level Numeric variable that specifies the significance level of adjusted p-value. A MATCH is >\code{sig.level}. Default: \code{0.1}.
#' @param reorder Boolean variable indicating if to reorder the columns so that matches are aligned along the diagonal. Default: \code{TRUE}.
#' @param return.value Boolean variable indicating if to return the plotted values. Default: \code{FALSE}.
#' @param main Plot title.
#' @param filename,width,height,units Plotting parameters passed to \code{\link[ggplot2]{ggsave}}.
#'
#' @return If \code{return.value = TRUE}, a matrix of \code{plotted.values}, and \code{pmat.adj} and \code{cell2cluster.adj} after adjustment.
#'
#' @export

plot_FRmatch_cell2cluster <- function(rst.cell2cluster, type="match", p.adj.method="BH", sig.level=0.1,
                                      reorder=TRUE, return.value=FALSE,
                                      main=NULL, filename=NA, width=NULL, height=NULL){

  pmat <- rst.cell2cluster$pmat #query cells in rows, ref cluster in columns
  pmat.adj <- apply(pmat, 1, function(x) p.adjust(x, method=p.adj.method)) %>% t() #row-wise p-value adjustment, i.e. per query cell
  ## update the scores with adjusted p-values, and determine matches and unassigned
  df <- rst.cell2cluster$cell2cluster %>% mutate(score=apply(pmat.adj,1,max)) %>%
    mutate(match=ifelse(score<sig.level, "unassigned", match))

  ## table of matches
  tab.match <- table(df$match, df$query.cluster) #query subclass in columns, prediction in rows
  if(!"unassigned" %in% rownames(tab.match)){
    tab.match <- rbind(tab.match,"unassigned"=0) #add the "unassigned" row if all matched
    tab.match <- as.table(tab.match) #still make it as a table class
  }
  tab.match.prop <- sweep(tab.match,2,colSums(tab.match),"/") #column sums should be 1

  ## order reference cluster orders in rows
  clusterNames.ref <- colnames(rst.cell2cluster$pmat) #reference cluster names IN ORDER
  oo <- match(c(clusterNames.ref,"unassigned"), rownames(tab.match.prop))
  oo.names <- rownames(tab.match.prop)[oo] %>% na.omit() #some ref clusters may not have matched query cells
  tab.match.prop <- tab.match.prop[oo.names,]

  ## reorder columns for query
  if(reorder) tab.match.prop <- reorder(tab.match.prop)
  oo.query.clusters <- colnames(tab.match.prop)

  ## plot: matching proportions bubble plot
  if(type=="match"){
    long.tab.match.prop <- tab.match.prop %>% as.data.frame() %>%
      select(query.cluster=Var2, match=Var1, Prop=Freq) %>%
      mutate(match=factor(match, levels = rev(c(clusterNames.ref, "unassigned"))))
    ## plot
    if(is.null(main)) main <- "FR-Match cell-to-cluster"
    g <- ggplot(long.tab.match.prop, aes(x=query.cluster, y=match, size=Prop, fill=Prop)) +
      geom_point(alpha = 0.7, shape = 21, color = "black") +
      scale_size_continuous(range = c(0, 10), limits = c(0, 1)) +
      scale_fill_viridis(option="D", guide = "legend", limits = c(0, 1)) +
      scale_y_discrete(drop=FALSE) + #show all ref clusters even if no match
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(main)
    ## save plot or plot on device
    if(is.null(width)) width <- ncol(tab.match.prop) + max(nchar(colnames(tab.match.prop)))/5 + 3
    if(is.null(height)) height <- length(clusterNames.ref) + max(nchar(clusterNames.ref))/5 + 1 #show all ref clusters even if no match
    if(!is.na(filename)){
      ggsave(filename, g, width = width, height = height, units = "cm")
    }
    else plot(g)
    ## output
    if(return.value){
      return(list("plotted.values"=tab.match.prop, "cell2cluster.adj"=df))
    }
  }

  ## plot: adjusted p-values boxplot
  if(type=="padj"){
    df.padj <- df %>% select(query.cluster, score) %>%
      mutate(query.cluster = fct_relevel(query.cluster, oo.query.clusters))
    g <- ggplot(df.padj, aes(x=query.cluster, y=score)) +
      geom_boxplot(outlier.size=.5) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      geom_hline(linetype = "dashed", yintercept = sig.level, color = "red") +
      scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
      xlab("Query cluster") + ylab("Adjusted p-value")
    ## save plot or plot on device
    if(!is.na(filename)) ggsave(filename, g, width=length(oo.query.clusters)*.2, height=5)
    else plot(g)
    ## output
    if(return.value) return(list("pmat.adj"=pmat.adj)) #plotted are the max padj for each query cell!!!
  }

  ## plot: force to match all
  if(type=="match.all"){
    df <- rst.cell2cluster$cell2cluster
    ## table of matches
    tab.match <- table(df$match, df$query.cluster) #query subclass in columns, prediction in rows
    tab.match.prop <- sweep(tab.match, 2, colSums(tab.match),"/") #column sums should be 1

    ## order reference cluster orders in rows
    clusterNames.ref <- colnames(rst.cell2cluster$pmat) #reference cluster names IN ORDER
    oo <- match(c(clusterNames.ref), rownames(tab.match.prop))
    oo.names <- rownames(tab.match.prop)[oo] %>% na.omit() #some ref clusters may not have matched query cells
    tab.match.prop <- tab.match.prop[oo.names,]

    ## reorder columns for query
    if(reorder) tab.match.prop <- reorder(tab.match.prop)

    ## plot
    long.tab.match.prop <- tab.match.prop %>% as.data.frame() %>%
      select(query.cluster=Var2, match=Var1, Prop=Freq) %>%
      mutate(match=factor(match, levels = rev(c(clusterNames.ref))))
    ## plot
    if(is.null(main)) main <- "FR-Match cell-to-cluster (match all)"
    g <- ggplot(long.tab.match.prop, aes(x=query.cluster, y=match, size=Prop, fill=Prop)) +
      geom_point(alpha = 0.7, shape = 21, color = "black") +
      scale_size_continuous(range = c(0, 10), limits = c(0, 1)) +
      scale_fill_viridis(option="D", guide = "legend", limits = c(0, 1)) +
      scale_y_discrete(drop=FALSE) + #show all ref clusters even if no match
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      ggtitle(main)
    ## save plot or plot on device
    if(is.null(width)) width <- ncol(tab.match.prop) + max(nchar(colnames(tab.match.prop)))/5 + 3
    if(is.null(height)) height <- length(clusterNames.ref) + max(nchar(clusterNames.ref))/5 + 1 #show all ref clusters even if no match
    if(!is.na(filename)){
      ggsave(filename, g, width = width, height = height, units = "cm")
    }
    else plot(g)
    ## output
    if(return.value){
      return(list("plotted.values"=tab.match.prop))
    }
  }
}
