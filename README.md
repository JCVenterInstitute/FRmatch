# FR-Match: A cluter-to-cluster cell type matching method for single cell RNA-seq experiments

## Getting Started

Installation:

```R
install.packages("devtools")
devtools::install_github("JCVenterInstitute/FRmatch", build_vignettes = TRUE)
```

To start with a demo Shiny App:

```R
FRmatch::runShiny()
```

To start with vignette:

```R
library(FRmatch)
browseVignettes("FRmatch")
```

### Prerequisites

* R, Shiny
* Data class: [Bioconductor-SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html)
* Feature selection: [JCVenterInstitute/NSForest](https://github.com/JCVenterInstitute/NSForest) (suggested)

### Description

Recently, the emergence of single cell RNA sequencing (scRNAseq) is providing large amounts of single cell transcriptomics data for the unbiased quantifications of cellular heterogeneity. Though scRNAseq data have been successfully generated by many labs, less attention has been paid to how knowledge derived from these data can be integrated across studies and leveraged by the whole single cell community.  In this R package, we provide a user-friendly scRNAseq integration tool that uses statistical methods to map new/query cell cluster data to the reference cell clusters.

Our method, FR-Match, is a novel application of the Friedman-Rafsky (FR) test, a non-parametric statistical test for multivariate data comparison in the context of single cell clustering results. We tailor the classical testing procedure for scRNAseq experiment data under the null hypothesis that there is no distributional difference in the two comparing clusters (i.e. a match) and the alternative hypothesis that the distributions of the two comparing clusters are different (i.e. a non-match) in the high-dimensional data space defined by selected gene features. Our procedure takes clustered gene expression matrices of query and reference experiments, and returns the FR statistic with adjusted p-value as evidence that the pair of comparing cell clusters is matched or not.

Manuscript submitted.

![](vignettes/FRmatch-scheme-v2.png)

## Versioning

This is version 0.0009 (initial release).
FRmatch is undergoing active development. Latest version and version control are managed in GitHub.

## Authors

* Yun (Renee) Zhang zhangy@jcvi.org
* Brian Aevermann baeverma@jcvi.org
* Richard Scheuermann RScheuermann@jcvi.org

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.

## Acknowledgments

* Allen Institute for Brain Science
* Chan Zuckerberg Initiative
