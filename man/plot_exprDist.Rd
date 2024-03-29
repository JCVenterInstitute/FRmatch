% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_exprDist.R
\name{plot_exprDist}
\alias{plot_exprDist}
\title{Gene expression data distribution plot}
\usage{
plot_exprDist(
  sce.E1,
  sce.E2,
  markers = NULL,
  name.E1 = "E1",
  name.E2 = "E2",
  breaks = 20,
  xlim = c(0, 10),
  ylim = c(0, 2),
  filename = NA,
  width = 10,
  height = 5,
  ...
)
}
\arguments{
\item{sce.E1, sce.E2}{Data objects, namely E1 and E2.}

\item{name.E1, name.E2}{Customized names for E1 and E2. Default: \code{"E1"} and \code{"E2"}, respectively.}

\item{breaks, xlim, ylim}{Plotting parameters passed to histogram plot.}

\item{filename}{File name if to save the plot. Default: \code{NA}, not to save the plot.}

\item{width, height}{Width and height for saved plot.}

\item{...}{Additional arguments passed to \code{\link[graphics]{hist}}}
}
\description{
This function plots the expression data distributions of the two single cell datasets (e.g. query and reference) to be compared.
}
