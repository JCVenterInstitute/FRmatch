## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  catch = TRUE,
  collapse = TRUE,
  comment = "#>"
)

## ----install, eval=FALSE------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("JCVenterInstitute/FRmatch")

## ----package, message=FALSE---------------------------------------------------
library(FRmatch)
library(SingleCellExperiment)
library(dplyr)
library(tibble)
library(data.table)

## ----data-example, warning=FALSE, message=FALSE-------------------------------
data(sce.example)
sce.example

## ----data-layer1--------------------------------------------------------------
## rename the data object
sce.layer1 <- sce.example
## cell type clusters and cluster sizes
knitr::kable(table(colData(sce.layer1)$cluster_membership), col.names=c("Cluster", "Size"))

## ----data-read----------------------------------------------------------------
## read in pieces of input data - this may take a few minutes
cell_by_gene_expression <- fread("cell_by_gene_expression.csv")
cell_cluster_labels <- fread("cell_cluster_labels.csv")
NSForest_marker_genes <- fread("NSForest_marker_genes.csv")
NSForest_fscores <- fread("NSForest_fscores.csv")
MTG_taxonomy <- fread("MTG_taxonomy.csv")$x #need to be a vector

## unique markers
unique_markers <- unique(NSForest_marker_genes$markerGene)

## ----make-data-object---------------------------------------------------------
sce.MTG <- make_data_object(dat = cell_by_gene_expression,
                            tab = cell_cluster_labels,
                            markers = unique_markers,
                            ## below are optional input data
                            cluster_marker_info = NSForest_marker_genes,
                            f_score = NSForest_fscores,
                            cluster_order = MTG_taxonomy)

## ----data-MTG-----------------------------------------------------------------
sce.MTG

## ----barcode-plot, fig.height = 5, fig.width = 7------------------------------
plot_cluster_by_markers(sce.example, cluster.name = "e1_e299_SLC17A7_L5b_Cdh13", name.self = "Layer1_")
plot_cluster_by_markers(sce.example, cluster.name = "g1_g48_GLI3_Astro_Gja1", name.self = "Layer1_")
plot_cluster_by_markers(sce.example, cluster.name = "i1_i90_COL5A2_Ndnf_Car4", name.self = "Layer1_")

## ----plot-clusterSize, fig.height = 9, fig.width = 9--------------------------
plot_clusterSize(sce.layer1, sce.MTG, name.E1 = "Layer1", name.E2 = "MTG")

## ----FRmatch_layer1toMTG, warning=FALSE, message=FALSE------------------------
rst.layer1toMTG <- FRmatch(sce.query = sce.layer1, sce.ref = sce.MTG, subsamp.size = 10)

## ----FRmatch_MTGtolayer1, warning=FALSE, message=FALSE------------------------
rst.MTGtolayer1 <- FRmatch(sce.query = sce.MTG, sce.ref = sce.layer1, subsamp.size = 10)

## ----FRmatch-bi-plot, fig.height = 13.5, fig.width = 7------------------------
plot_bi_FRmatch(rst.layer1toMTG, rst.MTGtolayer1, name.E1="Layer1_", name.E2="MTG_")

## ----FRmatch-plot-matches, fig.height = 13.5, fig.width = 7-------------------
## to visualize the one-directional matches
plot_FRmatch(rst.layer1toMTG)

## ----FRmatch-plot-padj, fig.height = 4, fig.width = 7-------------------------
## to visualize the adjusted p-values for each query cluster
plot_FRmatch(rst.layer1toMTG, type = "padj")

## ----FRmatch-plot-matches-2, fig.height = 5, fig.width = 13.5-----------------
## to visualize the one-directional matches
plot_FRmatch(rst.MTGtolayer1)

## ----FRmatch-plot-padj-2, fig.height = 4, fig.width = 8-----------------------
## to visualize the adjusted p-values for each query cluster
plot_FRmatch(rst.MTGtolayer1, type = "padj")

## ----FR-test, fig.height = 7, fig.width = 7-----------------------------------
# simulate some synthetic data from the same distribution
samp1 <- matrix(rnorm(1000), nrow = 50) #a 50-by-20 matrix: 50 dimensional, 20 data points
samp2 <- matrix(rnorm(1000), nrow = 50) #a 50-by-20 matrix: 50 dimensional, 20 data points
# FR test with MST plot
FRtest(samp1, samp2, plot.MST = TRUE, main = "Minimum spanning tree")

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

