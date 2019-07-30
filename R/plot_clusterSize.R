
#' Cluster size plot
#'
#' A function to plot cluster sizes for two datasets (e.g. query and reference).
#'
#' @param sce.E1,sce.E2 Two \code{FRmatch} input data objects, namely \code{E1} and \code{E2}. See example in \code{\link[FRmatch]{sce.example}}.
#' @param name.E1,name.E2 Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.
#' @param col.E1,col.E2 Customized colors for E1 and E2. Default: \code{"#F0E442"} and \code{"#56B4E9"}, respectively.
#' @param filename File name to save the plot. Default: \code{NA}, not to save.
#'
#' @export

plot_clusterSize <- function(sce.E1, sce.E2, name.E1="E1", name.E2="E2",
                             # filter.size=10,
                             col.E1="#F0E442", col.E2="#56B4E9", filename=NA){
  ## cluster sizes
  tab.E1 <- table(colData(sce.E1)$cluster_membership)
  tab.E2 <- table(colData(sce.E2)$cluster_membership)

  ## plot
  df.E1 <- tibble::tibble(cluster=names(tab.E1), size=tab.E1)
  g.E1 <- ggplot2::ggplot(data=df.E1, aes(x=cluster, y=size)) +
    ggplot2::geom_bar(stat="identity", position=position_dodge(), fill=col.E1) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggplot2::geom_text(aes(label=size), vjust=-0.3, size=2.5, position = position_dodge(0.9)) +
    # ggplot2::geom_hline(yintercept=filter.size, linetype="dashed", color="red") +
    ggplot2::ggtitle(name.E1)

  df.E2 <- tibble::tibble(cluster=names(tab.E2), size=tab.E2)
  g.E2 <- ggplot2::ggplot(data=df.E2, aes(x=cluster, y=size)) +
    ggplot2::geom_bar(stat="identity", position=position_dodge(), fill=col.E2) +
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggplot2::geom_text(aes(label=size), vjust=-0.3, size=2.5, position = position_dodge(0.9)) +
    # ggplot2::geom_hline(yintercept=filter.size, linetype="dashed", color = "red") +
    ggplot2::ggtitle(name.E2)

  g <- gridExtra::grid.arrange(g.E1, g.E2)
  plot(g)
  if(!is.na(filename)) ggplot2::ggsave(filename, g, width=max(length(tab.E1),length(tab.E2))*.2, height=10)
}



