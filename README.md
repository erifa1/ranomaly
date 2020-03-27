
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rANOMALY

<!-- badges: start -->

<!-- [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://www.tidyverse.org/lifecycle/#experimental) -->

<!-- badges: end -->

rANOMALY is the R Package of the first version of
[ANOMALY](https://forgemia.inra.fr/umrf/anomaly). [Here the
poster](https://prodinra.inra.fr/ft?id=%7BE7948F41-44A8-4F76-898D-92DDE09E40B9%7D&original=true)
presenting this workflow.

## Installation

You can install the development version of rANOMALY from this repository
with:

### Linux (highly recommended)

``` r
install.packages("devtools")
install.packages("remotes")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("biomformat")library(biomformat)

remotes::install_github("mahendra-mariadassou/phyloseq-extended")
devtools::install_git("https://forgemia.inra.fr/umrf/ranomaly")


# if install fails due to dependencies, retry: 
# remotes::install_github("mahendra-mariadassou/phyloseq-extended")
# devtools::install_git("https://forgemia.inra.fr/umrf/ranomaly")
```

### Windows

Require [Rtools](https://cran.r-project.org/bin/windows/Rtools/),
[git](https://git-scm.com/download/win) and run same commands as Linux
installation.
