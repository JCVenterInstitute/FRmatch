
#' Cluster size plot
#'
#' This is an auxilirary function that plots cluster sizes for two single cell clustered datasets (e.g. query and reference).
#'
#' @param sce.E1,sce.E2 Two \code{FRmatch} input data objects, namely \code{E1} and \code{E2}. See example in \code{\link[FRmatch]{sce.example}}.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param col.E1,col.E2 Customized colors for E1 and E2. Default: \code{"#F0E442"} and \code{"#56B4E9"}, respectively.
#' @param filename File name to save the plot. Default: \code{NA}, not to save.
#'
#' @export

plot_clusterSize <- function(sce.E1, sce.E2, name.E1="E1", name.E2="E2",
                             # filter.size=10,
                             col.E1="#F0E442", col.E2="#56B4E9",
                             filename=NA, width=NULL, height=10){
  ## cluster sizes
  tab.E1 <- table(colData(sce.E1)$cluster_membership)
  tab.E2 <- table(colData(sce.E2)$cluster_membership)

  ## plot
  df.E1 <- tibble::tibble(Cluster=names(tab.E1), Size=tab.E1)
  g.E1 <- ggplot2::ggplot(data=df.E1, ggplot2::aes(x=Cluster, y=Size)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge(), fill=col.E1) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::geom_text(ggplot2::aes(label=Size), vjust=-0.3, size=2.5, position = ggplot2::position_dodge(0.9)) +
    # ggplot2::geom_hline(yintercept=filter.size, linetype="dashed", color="red") +
    ggplot2::ggtitle(paste0(name.E1, " (", length(tab.E1)," clusters)"))

  df.E2 <- tibble::tibble(Cluster=names(tab.E2), Size=tab.E2)
  g.E2 <- ggplot2::ggplot(data=df.E2, ggplot2::aes(x=Cluster, y=Size)) +
    ggplot2::geom_bar(stat="identity", position=ggplot2::position_dodge(), fill=col.E2) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::geom_text(ggplot2::aes(label=Size), vjust=-0.3, size=2.5, position = ggplot2::position_dodge(0.9)) +
    # ggplot2::geom_hline(yintercept=filter.size, linetype="dashed", color = "red") +
    ggplot2::ggtitle(paste0(name.E2, " (", length(tab.E2)," clusters)"))

  g <- gridExtra::grid.arrange(g.E1, g.E2)
  plot(g)

  ## save plot
  if(!is.na(filename)){
    if(is.null(width)) width <- max(length(tab.E1), length(tab.E2))*.2
    ggplot2::ggsave(filename, g, width=width, height=height)
  }
}



