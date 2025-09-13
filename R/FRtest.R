
#' Friedman-Rafsky (FR) test
#'
#' FR test is a multivariate generalization of nonparametric two-sample test. This function is an implementation
#' with customized options, including a visualization of the minimum spanning tree (MST).
#'
#' @param samp1 Numeric matrix or data frame for Sample 1. Rows are multivariate dimensions, and columns are samples. E.g. gene by cell.
#' @param samp2 Numeric matrix or data frame for Sample 2.
#' @param use.cosine An option if to use cosine distance. Logical variable. By default (\code{FALSE}),
#' Euclidean distance is used.
#' @param binary An option if to use binary values. Logical variable. Default: \code{FALSE}. If \code{TRUE},
#' use \code{binary.cutoff} to dichotomize \code{samp1} and \code{samp2}.
#' @param binary.cutoff Numeric value for binary cutoff. Binary value = 1 if greater than \code{binary.cutoff}, 0 otherwise. Default: \code{2}.
#' @param plot.MST Boolean variable indicating if to plot the minimum spanning tree (MST). Default: \code{FALSE}.
#' @param col Character vector of length two for customized colors of the nodes in MST. Default: \code{c("#F0E442", "#56B4E9")}.
#' @param label.names Character vector of length two for customized names of the two samples. Default: \code{c("Sample 1","Sample 2")}.
#' @param vertex.size,edge.width,... Additional plotting parammeters passed to \code{\link[igraph]{plot.igraph}}. Default: \code{vertex.size=5, edge.width=1}.
#'
#' @return Test statistics and p-values.
#' \item{runs}{Total number of subtrees.}
#' \item{runs.samp1}{Number of subtrees of Sample 1.}
#' \item{runs.samp2}{Number of subtrees of Sample 2.}
#' \item{stat}{The standardized FR statistic.}
#' \item{p.value}{P-value of the FR test.}
#'
#' @examples
#' \dontrun{
#' samp1 <- matrix(rnorm(100),nrow=5)
#' samp2 <- matrix(rnorm(100),nrow=5)
#' FRtest(samp1, samp2)
#' FRtest(samp1, samp2, use.cosine=TRUE)
#' FRtest(samp1, samp2, plot.MST=TRUE, main="Minimum spanning tree plot")
#' FRtest(samp1, samp2, binary=TRUE, binary.cutoff=1)
#' }
#' @export

FRtest <- function(samp1, samp2,
                    ## two options: 1. cosine distance, 2. use binary values
                    use.cosine=FALSE, binary=FALSE, binary.cutoff=2,
                    ## plot minimum spanning tree
                    plot.MST=FALSE, col=c("#F0E442", "#56B4E9"), label.names=c("Sample 1","Sample 2"),
                    vertex.size=5, edge.width=1,
                    ...)
{
  ## data input matrices: rows = multivariate dimensions, columns = samples
  xx <- as.matrix(samp1)
  yy <- as.matrix(samp2)

  ## get sample sizes
  m <- max(ncol(xx), 1)
  n <- max(ncol(yy), 1)
  N <- m+n

  ##--- OPTION 2: BINARY VALUES ---##
  if(binary){
    xx <- matrix(as.numeric(xx>binary.cutoff), nrow=nrow(xx))
    yy <- matrix(as.numeric(yy>binary.cutoff), nrow=nrow(yy))
  }
  ##---##

  ## column-wise combined data matrix
  dat <- cbind(xx, yy)

  ##--- OPTION 1: CONSINE DISTANCE ---##
  if(use.cosine){
    cosmat <- lsa::cosine(dat)
    cosmat[is.na(cosmat)] <- -1 #if one point is at the origin, set the -1,0,1??? cosine value
    distobj <- as.dist(1-cosmat) #large cosine, small distance; and vise versa
  }
  ##---##

  ## euclidean distance and MST
  else distobj <- dist(t(dat), method="euclidean", diag=T, upper=T)
  myMST <- ade4::neig2mat(ade4::mstree(distobj))


  ## count total runs (i.e. total subgraphs)
  bottomleft <- myMST[(m+1):N,1:m]
  runs <- sum(bottomleft)+1

  ## count subgraphs for each sample
  bottomright <- myMST[(m+1):N,(m+1):N]
  g.samp2 <- igraph::graph.adjacency(bottomright, mode="upper", weighted=NULL, diag=FALSE)
  runs.samp2 <- igraph::components(g.samp2)$no
  topleft <- myMST[1:m,1:m]
  g.samp1 <- igraph::graph.adjacency(topleft, mode="upper", weighted=NULL, diag=FALSE)
  runs.samp1 <- igraph::components(g.samp1)$no
  ## check
  # runs==sum(runs.samp1, runs.samp2)

  ## calculate common nodes
  xsum <- colSums(myMST)
  C <- sum(xsum*(xsum-1))/2
  ## calculate mean and variance
  mu <- 2*m*n/N + 1
  sigma.sq <- (2*m*n/(N*(N-1)))*((2*m*n-N)/N+(C-N+2)*(N*(N-1)-4*m*n+2)/((N-2)*(N-3)))
  ## the standardized FR-stat and p-value
  stat <- (runs-mu)/sqrt(sigma.sq)
  p.value <- pnorm(stat)

  ## plot
  if(plot.MST){
    # adjmat <- as.matrix(distobj)*myMST
    # g <- igraph::graph.adjacency(adjmat, mode="upper", weighted=TRUE, diag=FALSE)
    g <- igraph::graph.adjacency(myMST, mode="upper", weighted=NULL, diag=FALSE)
    colors <- rep(col, c(m,n))
    plot(g, vertex.size=vertex.size, edge.width=edge.width, vertex.label=NA, vertex.color=colors, frame=TRUE, ...)
    legend("bottom", paste0(label.names, " (",c(m,n),")"),
                            fill=col, bty ="n", xpd=TRUE, inset=c(0, -.18))
  }

  ## output
  return(c("runs"=runs, "runs.samp1"=runs.samp1, "runs.samp2"=runs.samp2, "stat"=stat, "p.value"=p.value))
}



