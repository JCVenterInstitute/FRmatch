# #' Friedman-Rafsky (FR) test with subsampling
# #'
# #' FR test with subsampling, which is provided for comparing two samples with one or both large sample size(s).
# #'
# #' @param samp1 Numeric matrix or data frame for Sample 1.
# #' @param samp2 Numeric matrix or data frame for Sample 2.
# #' @param subsamp.size,subsamp.iter Cluster size and number of iterations for subsampling.
# #' @param ... Additional parammeters passed to \code{\link[FRmatch]{FR.test}}.
# #'
# #' @return Median p-value from the subsampling itertations and its corresponding test statistics will be returned.
# #' \item{runs}{Total number of subtrees.}
# #' \item{runs.samp1}{Number of subtrees of Sample 1.}
# #' \item{runs.samp2}{Number of subtrees of Sample 2.}
# #' \item{stat}{The standardized FR statistic.}
# #' \item{p.value}{P-value of the FR test.}
# #'
# #' @examples
# #' \dontrun{
# #' samp1 <- matrix(rnorm(1000),nrow=5)
# #' samp2 <- matrix(rnorm(1000),nrow=5)
# #' FR.test.subsamp.each(samp1, samp2)
# #' }
#'
#' @importFrom dplyr %>%

FR.test.subsamp.each <- function(samp1, samp2, subsamp.size, subsamp.iter, ...){
  out.all <- NULL
  for(b in 1:subsamp.iter){
    ## data input matrices: rows = multivariate dimensions, columns = samples
    xx <- as.matrix(samp1)
    yy <- as.matrix(samp2)
    ## get original sample sizes
    m <- max(ncol(xx), 1)
    n <- max(ncol(yy), 1)
    ## subsampling sizes
    mm <- sample(1:m, min(subsamp.size,m), replace=FALSE)
    nn <- sample(1:n, min(subsamp.size,n), replace=FALSE)
    xx.B <- samp1[,mm]
    yy.B <- samp2[,nn]
    ## FR test on subsampled samples
    out.B <- FR.test(xx.B, yy.B, ...)
    out.all <- rbind(out.all, out.B)
  }
  out.all.sort <- out.all %>% as.data.frame() %>% dplyr::arrange(p.value)
  output <- out.all.sort[round(subsamp.iter/2),]
  output <- as.numeric(output)
  names(output) <- names(out.all.sort)
  return(output)
}

FR.test.subsamp <- function(samp1, samp2, subsamp.size, subsamp.iter, ...){
  ## first pass of FR test to check if it is a trivial seperation, i.e. one sample only forms one subtree
  out0 <- FR.test(samp1, samp2, ...)
  if(out0["runs.samp1"]==1 | out0["runs.samp2"]==1) output <- out0
  ## if not the trivial case, do subsampling
  else output <- FR.test.subsamp.each(samp1, samp2, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, ...)
  ## output
  return(output)
}


