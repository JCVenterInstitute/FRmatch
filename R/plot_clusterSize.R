
#' Cluster size plot
#'
#' This function plots the cluster sizes of the two single cell datasets (e.g. query and reference) to be compared.
#'
#' @param sce.E1,sce.E2 Data objects, namely E1 and E2.
#' @param decreasing Boolean variable if the cluster sizes should be sorted in decreasing order. Default: \code{TRUE}.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param col.E1,col.E2 Customized colors for E1 and E2. Default: \code{"#F0E442"} and \code{"#56B4E9"}, respectively.
#' @param filename File name if to save the plot. Default: \code{NA}, not to save the plot.
#' @param width,height Width and height for saved plot.
#'
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#'
#' @export

plot_clusterSize <- function(sce.E1, sce.E2, decreasing=TRUE,
                             # filter.size=10,
                             name.E1="E1", name.E2="E2", col.E1="#F0E442", col.E2="#56B4E9",
                             filename=NA, width=NULL, height=10){
  ## cluster sizes
  tab.E1 <- table(colData(sce.E1)$cluster_membership)
  tab.E2 <- table(colData(sce.E2)$cluster_membership)
  if(decreasing){
    tab.E1 <- sort(tab.E1, decreasing=TRUE)
    tab.E2 <- sort(tab.E2, decreasing=TRUE)
  }

  ## plot
  df.E1 <- tibble(Cluster=names(tab.E1), Size=as.vector(tab.E1))
  g.E1 <- ggplot(data=df.E1, aes(x=Cluster, y=Size)) +
    geom_bar(stat="identity", position=position_dodge(), fill=col.E1) +
    scale_x_discrete(limits = names(tab.E1)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    geom_text(aes(label=Size), vjust=-0.3, size=2.5, position=position_dodge(0.9)) +
    # geom_hline(yintercept=filter.size, linetype="dashed", color="red") +
    ggtitle(paste0(name.E1, " (", length(tab.E1)," clusters, ",sum(tab.E1)," cells)"))

  df.E2 <- tibble(Cluster=names(tab.E2), Size=as.vector(tab.E2))
  g.E2 <- ggplot(data=df.E2, aes(x=Cluster, y=Size)) +
    geom_bar(stat="identity", position=position_dodge(), fill=col.E2) +
    scale_x_discrete(limits = names(tab.E2)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=90, hjust=1)) +
    geom_text(aes(label=Size), vjust=-0.3, size=2.5, position=position_dodge(0.9)) +
    # geom_hline(yintercept=filter.size, linetype="dashed", color="red") +
    ggtitle(paste0(name.E2, " (", length(tab.E2)," clusters, ",sum(tab.E2)," cells)"))

  g <- grid.arrange(g.E1, g.E2)
  plot(g)

  ## save plot
  if(!is.na(filename)){
    if(is.null(width)) width <- max(length(tab.E1), length(tab.E2))*.2
    ggsave(filename, g, width=width, height=height)
  }
}



