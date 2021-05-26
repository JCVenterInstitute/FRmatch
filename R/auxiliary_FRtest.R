
########################################################################################################
## FRtest_subsamp()
########################################################################################################

FRtest_subsamp <- function(samp1, samp2, subsamp.size, subsamp.iter, return.all=FALSE, ...){
  ## data input matrices: rows = multivariate dimensions, columns = samples
  xx <- as.matrix(samp1)
  yy <- as.matrix(samp2)
  ## get original sample sizes
  m <- max(ncol(xx), 1)
  n <- max(ncol(yy), 1)

  out.all <- NULL
  for(b in 1:subsamp.iter){
    ## subsampling sizes
    mm <- sample(1:m, min(subsamp.size,m), replace=FALSE)
    nn <- sample(1:n, min(subsamp.size,n), replace=FALSE)
    xx.B <- samp1[,mm]
    yy.B <- samp2[,nn]
    ## FR test on subsampled samples
    out.B <- FRtest(xx.B, yy.B, ...)
    out.all <- rbind(out.all, out.B)
  }
  out.all.sort <- out.all %>% as.data.frame() %>% arrange(p.value)
  if(return.all){
    return(out.all.sort)
  } else
  output <- out.all.sort[round(subsamp.iter/2),]
  output <- as.numeric(output)
  names(output) <- names(out.all.sort)
  return(output)
}

# FRtest_subsamp <- function(samp1, samp2, subsamp.size, subsamp.iter, return.all=FALSE, ...){
#   ## first pass of FR test to check if it is a trivial separation, i.e. one sample only forms one subtree
#   out0 <- FRtest(samp1, samp2, ...)
#   if(out0["runs.samp1"]==1 | out0["runs.samp2"]==1) output <- out0
#   ## if not the trivial case, do subsampling
#   else output <- FRtest_subsamp_each(samp1, samp2, subsamp.size=subsamp.size, subsamp.iter=subsamp.iter, return.all=return.all, ...)
#   ## output
#   return(output)
# }


########################################################################################################
## FRtest_cell2cluster()
########################################################################################################

FRtest_cell2cluster <- function(samp1, samp2, subsamp.size, subsamp.iter, ...){

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
    out.B <- FRtest(xx.B, yy.B, ...)
    out.cell2cluster[mm] %<>% pmax(out.B["p.value"], na.rm=TRUE) #pmax: parallel maxima, meaning replace the original value in out.cell2cluster if out.B["p.value"] is larger
  }
  return(out.cell2cluster)
}
