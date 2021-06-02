
# MicroVis

<!-- badges: start -->
<!-- badges: end -->

MicroVis is a package for flexible analysis of metagenomic data and generation of customizable, publication-ready figures

## Installation

Before installing MicroVis, make sure to install these dependencies that separately for full functionality:

Install [phyloseq](https://joey711.github.io/phyloseq/) and [DESeq2](http://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html):

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
```

Install [microbiomemarker]():
``` r
if (!require("devtools")) install.packages("devtools")

devtools::install_github("yiluheihei/microbiomeMarker") 
```

You can then install the released version of MicroVis from [GitHub](https://github.com/) with:

``` r
devtools::install_github("https://github.com/microresearcher/MicroVis.git")
```
