% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_bi_FRmatch.R
\name{plot_bi_FRmatch}
\alias{plot_bi_FRmatch}
\title{Plotting function for bi-directional FR-Match results}
\usage{
plot_bi_FRmatch(
  rst.FRmatch.E1toE2,
  rst.FRmatch.E2toE1,
  prefix = c("query.", "ref."),
  name.E1 = "E1.",
  name.E2 = "E2.",
  p.adj.method = "BY",
  sig.level = 0.05,
  reorder = TRUE,
  two.way.only = FALSE,
  return.value = FALSE,
  cellwidth = 10,
  cellheight = 10,
  main = NULL,
  filename = NA,
  ...
)
}
\arguments{
\item{rst.FRmatch.E1toE2, rst.FRmatch.E2toE1}{The \code{\link[FRmatch]{FRmatch}} outputs.}

\item{prefix}{Prefix names for query and reference clusters used in \code{\link[FRmatch]{FRmatch}}.}

\item{name.E1, name.E2}{Customized names with delimiter for E1 and E2 to be used in this figure. Default: \code{"E1."} and \code{"E2."}, respectively.}

\item{p.adj.method}{See \code{\link[FRmatch]{plot_FRmatch}}.}

\item{sig.level}{See \code{\link[FRmatch]{plot_FRmatch}}.}

\item{reorder}{See \code{\link[FRmatch]{plot_FRmatch}}.}

\item{two.way.only}{Boolean variable indicating if to plot two-way matches only. Default: \code{FALSE}.}

\item{return.value}{Boolean variable indicating if to return the plotted values. Default: \code{FALSE}.}

\item{cellwidth, cellheight, main, filename, ...}{Plotting parameters passed to \code{\link[pheatmap]{pheatmap}}.}
}
\value{
If \code{return.value = TRUE}, a matrix of two-way matching values 2 = two-way match, 1 = one-way match, and 0 = no match.
}
\description{
This function combines and plots two sets of reciprocal \code{\link[FRmatch]{FRmatch}} results (i.e. experiment 1 (E1) query
to experiment 2 (E2) reference with E2 markers, and E2 query to E1 reference with E1 markers).
}
\seealso{
\code{\link[FRmatch]{plot_FRmatch}}.
}
