




FR.test.B <- function(samp1, samp2, Biter=1000, seed=NULL, ...)
{
  ## data input matrices: rows = multivariate dimensions, columns = samples
  xx <- as.matrix(samp1)
  yy <- as.matrix(samp2)

  ## get sample sizes
  m <- max(ncol(xx), 1)
  n <- max(ncol(yy), 1)
  N <- m+n

  ## column-wise combined data matrix
  dat <- cbind(xx,yy)

  ## BOOTSTRAPED / PERMUTED? data <--- THIS IS WRONG!!!!!! JUST KEEP THE FORMAT AND CORRECT IT LATER.
  set.seed(seed)
  ind.B <- lapply(1:Biter, function(b) sample(1:N)) #permutate combined data points
  dat.B <- lapply(ind.B, function(ind.b) dat[,ind.b])

  ## FR test
  # out.B <- lapply(dat.B, function(dat.b) FR.test(dat.b[,1:m], dat.b[,1:m], plot.MST=FALSE, ...))
  out.B <- lapply(dat.B, function(dat.b) FR.test(dat.b[,1:m], dat.b[,1:m], plot.MST=FALSE))
  out.B <- do.call("rbind", out.B)
  out <- out.B %>% as.data.frame() %>% arrange(p.value) %>% slice(ceiling(Biter/2)) %>% unlist()

  ## return
  return(out)
}

# samp1 <- matrix(rnorm(1000),nrow=5)
# samp2 <- matrix(rnorm(100,1),nrow=5)
# FR.test(samp1, samp2)
# FR.test.B(samp1, samp2)
# FR.test.B(samp1, samp2, seed=1)
# FR.test.B(samp1, samp2, use.cosine = TRUE, seed=1)
# FR.test(samp1, samp2, col=c("#F0E442", "#0072B2"), label.names = c("qeury", "ref"), main="Minimum spanning tree plot")
# FR.test(samp1, samp2, binary = TRUE, binary.cutoff = 1, seed=1)



FR.test.boot <- function(xx, yy, N.samp=100, K=20, Biter=1000, use.cosine=FALSE)
{
  m <- max(dim(xx)[2],1)
  n <- max(dim(yy)[2],1)
  N <- sum(m,n)

  if(N>N.samp){ #<------is this a good criterion to kick off bootstrap? Is this the problem for FR-test unstability?
    ## BOOTSTRAP
    # N.samp <- 100
    # K <- 20
    pvals.all <- NULL
    for(k in 1:K){
      # cat(k,"\n")
      mm <- sample(1:m, max(floor(m/N*N.samp),1), replace=FALSE) #<------ WITHOUT REPLACEMENT!!! OR NOT???
      nn <- sample(1:n, max(floor(n/N*N.samp),1), replace=FALSE)
      xx.B <- xx[,mm]
      yy.B <- yy[,nn]
      out.org <- FR.test(xx.B, yy.B, plot.MST=FALSE, use.cosine=use.cosine)

      out.B <- NULL
      # Biter <- 1000
      for(B in 1:Biter){
        out <- FR.test.B(xx.B, yy.B, plot.MST=FALSE, use.cosine=use.cosine)
        out.B <- rbind(out.B, out)
      }
      # out.all.sort <- as.matrix(tbl_df(out.all) %>% arrange(runs))
      pval.k <- sum(out.B[,"runs"]<=out.org["runs"])/Biter
      pvals.all <- c(pvals.all, pval.k)
      p.out <- median(pvals.all)
    }
  }
  else p.out <- FR.test(xx,yy)["p.value"]
  return(p.out)
}


# library(GSAR)
# library(MASS)
# ngenes <- 20
# nsamples <- 80
# ## let the mean vector have zeros of length 20 for condition 1
# zero_vector <- array(0,c(1,ngenes))
# ## let the mean vector have 2s of length 20 for condition 2
# mu_vector <- array(1,c(1,ngenes))
# ## set the covariance matrix to be an identity matrix
# cov_mtrx <- diag(ngenes)
# gp1 <- mvrnorm((nsamples/2), zero_vector, cov_mtrx)
# gp2 <- mvrnorm((nsamples/2), mu_vector, cov_mtrx)
# ## combine the data of two conditions into one dataset
# gp <- rbind(gp1,gp2)
# dataset <- aperm(gp, c(2,1))
# ## first 20 samples belong to group 1
# ## second 20 samples belong to group 2
# result <- WWtest(object=dataset, group=c(rep(1,nsamples/2),rep(2,nsamples/2)))
# result
#
#
# FR.test(dataset[,1:(nsamples/2)], dataset[,(nsamples/2+1):nsamples])
