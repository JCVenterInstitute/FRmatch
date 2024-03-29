% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_data_object.R
\name{make_data_object}
\alias{make_data_object}
\title{Make data object}
\usage{
make_data_object(
  dat,
  tab,
  markers,
  cluster_marker_info = NULL,
  f_score = NULL,
  cluster_order = NULL
)
}
\arguments{
\item{dat}{Cell-by-gene expression data in a data frame or equivalent format. The first column is cell names, followed by columns of each gene.}

\item{tab}{Cluster membership of each cell in a data frame or equivalent format. The first column is cell names, and the second column is cluster labels.}

\item{markers}{A vector of marker genes.}

\item{cluster_marker_info}{Optionally, a data frame of maker genes for each cluster. \code{marker} is the unique set of marker genes for all clusters.}

\item{f_score}{Optionally, a data frame of F-beta scores for each cluster.}

\item{cluster_order}{Optionally, a vector of ordered cluster names, which may reflect the taxonomy of the cell type clusters.}
}
\value{
A data object of the \link[SingleCellExperiment]{SingleCellExperiment} class
}
\description{
This function creates the data object used in \code{\link[FRmatch]{FRmatch}}.
See an example of the data object in \code{\link[FRmatch]{sce.example}}.
}
