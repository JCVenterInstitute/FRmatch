
#' Plotting function for FR-Match cell2cluster results
#'
#' @export

plot_FRmatch_cell2cluster <- function(rst.cell2cluster, type="match.prop", p.adj.method="BH", sig.level=0.05,
                                      reorder=TRUE, return.value=FALSE,
                                      filename=NA, width=NULL, height=NULL){

  ## calculate adjusted p-values
  clusterNames.ref <- names(rst.cell2cluster)[-(1:4)] #reference cluster names IN ORDER
  pmat <- rst.cell2cluster[,clusterNames.ref]
  pmat.adj <- pmat %>% as_tibble() %>%
    mutate(across(everything(), ~ p.adjust(.x, method=p.adj.method)))
  rst.cell2cluster.padj <- rst.cell2cluster %>% select(query.cell, query.cluster) %>%
    mutate(match.cell2cluster=clusterNames.ref[max.col(pmat.adj)], rmax.cell2cluster=apply(pmat.adj,1,max)) %>%
    cbind(pmat.adj)

  ## determine matches and unassigned
  df.cell2cluster <- rst.cell2cluster.padj %>%
    select(query.cell, query.cluster, match.cell2cluster, rmax.cell2cluster) %>%
    mutate(match.cell2cluster=ifelse(rmax.cell2cluster<sig.level, "unassigned", match.cell2cluster))

  ### plot matching proportion of cells by query clusters ###
  if(type=="match.prop"){
    ## matching proportion bubble plot
    tab.match <- table(df.cell2cluster$match.cell2cluster, df.cell2cluster$query.cluster) #query subclass in columns, prediction in rows
    tab.match.prop <- sweep(tab.match,2,colSums(tab.match),"/") #column sums should be 1

    ## order reference cluster orders in rows
    oo <- match(c(clusterNames.ref,"unassigned"), rownames(tab.match.prop))
    oo.names <- rownames(tab.match.prop)[oo] %>% na.omit() #some ref clusters may not have matched query cells
    tab.match.prop <- tab.match.prop[oo.names,]

    ## reorder columns for query
    if(reorder) tab.match.prop <- reorder(tab.match.prop)

    ## for ggplot
    long.tab.match.prop <- tab.match.prop %>% as.data.frame() %>%
      select(query.cluster=Var2, match.cell2cluster=Var1, Prop=Freq) %>%
      mutate(match.cell2cluster=factor(match.cell2cluster, levels = rev(c(clusterNames.ref, "unassigned"))))

    ## plot
    g <- ggplot(long.tab.match.prop, aes(x=query.cluster, y=match.cell2cluster, size=Prop, fill=Prop)) +
      geom_point(alpha=0.7, shape=21, color="black") +
      scale_size_continuous(range = c(0, 10)) +
      scale_fill_viridis(option="D", guide = "legend") +
      scale_y_discrete(drop=FALSE) + #show all ref clusters even if no match
      theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    ## save plot or plot on device
    if(is.null(width)) width <- ncol(tab.match.prop)*.2+1
    if(is.null(height)) height <- nrow(tab.match.prop)*.2
    if(!is.na(filename)) ggsave(filename, g, width=width, height=height)
    else plot(g)

    ## output
    if(return.value){
      return("rst.cell2cluster.padj"=rst.cell2cluster.padj, "match.prop"=tab.match.prop)
    }
  }

  ### plot jaccard similarity between query and ref clusters ###

}
