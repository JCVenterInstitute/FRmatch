## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  catch = TRUE,
  collapse = TRUE,
  comment = "##"
)

## ---- install, eval=FALSE------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("JCVenterInstitute/FRmatch")

## ---- package, message=FALSE---------------------------------------------
library(FRmatch)

## ---- shinyapp, eval=FALSE-----------------------------------------------
#  runShiny()

## ---- data, warning=FALSE, message=FALSE---------------------------------
data(sce.example)
sce.example

## ---- cluster-size, , warning=FALSE, message=FALSE-----------------------
library(SingleCellExperiment)
knitr::kable(table(colData(sce.example)), col.names=c("Cluster", "Size"))

## ---- data-split, warning=FALSE, message=FALSE---------------------------
library(dplyr)
library(tibble)
set.seed(999)
## subsampling
all <- colData(sce.example) %>% as.data.frame() %>% rownames_to_column()
sam1 <- all %>% group_by(cluster_membership) %>% sample_frac(.5)
sam2 <- dplyr::setdiff(all, sam1)

sce.sam1 <- sce.example[,sam1$rowname] #query
sce.sam2 <- sce.example[,sam2$rowname] #reference

## ---- FRmatch, warning=FALSE, message=FALSE------------------------------
rst <- FRmatch(sce.query = sce.sam1, sce.ref = sce.sam2)

## ---- FRmatch-plot, fig.height=5, fig.width=6----------------------------
plot_FRmatch(rst)
plot_FRmatch(rst, type="padj")

## ---- dropout, fig.height=8, fig.width=7---------------------------------
plot_nonzero(sce.example, return.value=FALSE, return.plot=TRUE)

## ---- FR-test, fig.height=5, fig.width=5---------------------------------
samp1 <- matrix(rnorm(50),nrow=5)
samp2 <- matrix(rnorm(100),nrow=5)
FR.test(samp1, samp2, plot.MST=TRUE, main="Minimum spanning tree plot")

## ---- sessionInfo--------------------------------------------------------
sessionInfo()

