
FR.test.cell2cluster <- function(samp1, samp2, subsamp.size, subsamp.iter, ...){

  ## data input matrices: rows = multivariate dimensions, columns = samples
  xx <- as.matrix(samp1)
  yy <- as.matrix(samp2)
  ## get original sample sizes
  m <- max(ncol(xx), 1)
  n <- max(ncol(yy), 1)

  out.cell2cluster <- rep(NA, m)
  names(out.cell2cluster) <- colnames(xx)
  for(b in 1:subsamp.iter){
    ## subsampling sizes
    mm <- sample(1:m, min(subsamp.size,m), replace=FALSE)
    nn <- sample(1:n, min(subsamp.size,n), replace=FALSE)
    xx.B <- samp1[,mm]
    yy.B <- samp2[,nn]
    ## FR test on subsampled samples
    # out.B <- FR.test(xx.B, yy.B, ...)
    out.B <- FR.test(xx.B, yy.B)
    out.cell2cluster[mm] %<>% pmax(out.B["p.value"], na.rm=TRUE)
  }
  return(out.cell2cluster)
}
