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
## other useful packages for this vignette
library(SingleCellExperiment)
library(dplyr)
library(tibble)

## ---- shinyapp, eval=FALSE-----------------------------------------------
#  runShiny()

## ---- data, warning=FALSE, message=FALSE---------------------------------
data(sce.example)
sce.example

## ---- cluster-size, , warning=FALSE, message=FALSE-----------------------
knitr::kable(table(colData(sce.example)), col.names=c("Cluster", "Size"), row.names=1:15)

## ---- data-split, warning=FALSE, message=FALSE---------------------------
set.seed(999)
all <- colData(sce.example) %>% as.data.frame() %>% rownames_to_column()
sam1 <- all %>% group_by(cluster_membership) %>% sample_frac(.5)
sam2 <- dplyr::setdiff(all, sam1)

sce.sam1 <- sce.example[,sam1$rowname] #query
sce.sam2 <- sce.example[,sam2$rowname] #reference

## ---- clusterSize, fig.height = 9, fig.width = 7-------------------------
plot_clusterSize(sce.sam1, sce.sam2)

## ---- FRmatch, warning=FALSE, message=FALSE------------------------------
rst <- FRmatch(sce.query = sce.sam1, sce.ref = sce.sam2)

## ---- barcoding-plot-----------------------------------------------------
plot_cluster_by_markers(sce.sam1, cluster.name = "i1_i90_COL5A2_Ndnf_Car4")

## ---- FRmatch-plot-------------------------------------------------------
plot_FRmatch(rst)
plot_FRmatch(rst, type="padj")

## ---- FR-test------------------------------------------------------------
# simulate some simple data
samp1 <- matrix(rnorm(200),nrow=5) #a 5-by-40 matrix
samp2 <- matrix(rnorm(100),nrow=5) #a 5-by-20 matrix
# FR test with MST plot
FR.test(samp1, samp2, plot.MST=TRUE)

## ---- sessionInfo--------------------------------------------------------
sessionInfo()

