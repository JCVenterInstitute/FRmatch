% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_MST.R
\name{plot_MST}
\alias{plot_MST}
\title{Plot minimum spanning tree (MST)}
\usage{
plot_MST(sce.query, sce.ref, query.cluster, ref.cluster, nsamp = 30, ...)
}
\arguments{
\item{sce.query, sce.ref}{Query and reference data objects.}

\item{query.cluster, ref.cluster}{Query and reference cluster names to plot.}

\item{nsamp}{Number of randomly selected cells to plot for large cluster. Default: 30.}

\item{...}{Additional arguments passed to \code{\link[FRmatch]{FRtest}}.}
}
\value{
MST plot and FR-test result in console.
}
\description{
This function is a wrapper function for plotting MST of two interested clusters.
}
