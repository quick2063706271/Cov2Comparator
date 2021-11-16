
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Cov2Comparator

<!-- badges: start -->
<!-- badges: end -->

The goal of Cov2Comparator is to analyze SARS covid-2 genome across
different geographic regions including performing multiple sequence
alignments, phylogenetic tree and mutation analysis.

## Installation

You can install the development version of Cov2Comparator from
[GitHub](https://github.com/) with:

``` r
require("devtools")
devtools::install_github("quick2063706271/Cov2Comparator", build_vignettes = TRUE)
library("Cov2Comparator")
```

To run the shinyApp: Under construction

## Overview

``` r
ls("package:Cov2Comparator")
data(package = "Cov2Comparator") 
```

### Vignettes

``` r
browseVignettes("Cov2Comparator")
```

## Reference

R Core Team (2021). R: A language and environment for statistical
computing. R Foundation for Statistical Computing, Vienna, Austria. URL
<https://www.R-project.org/>.

Charif, D. and Lobry, J.R. (2007). SeqinR 1.0-2: a contributed package
to the R project for statistical computing devoted to biological
sequences retrieval and analysis. Structural approaches to sequence
evolution: Molecules, networks, populations, series Biological and
Medical Physics, Biomedical Engineering, 207-232. Springer Verlag, New
York.

Ulrich Bodenhofer, Enrico Bonatesta, Christoph Horejs-Kainrath, & Sepp
Hochreiter (2015). msa: an R package for multiple sequence alignment.
Bioinformatics, 31(24), 3997–3999.

Paradis E. & Schliep K. (2019). ape 5.0: an environment for modern
phylogenetics and evolutionary analyses in R. Bioinformatics 35:
526-528.

Raivo Kolde (2019). pheatmap: Pretty Heatmaps. R package version 1.0.12.
<https://CRAN.R-project.org/package=pheatmap>

## Acknowledgemments

This package was developed as part of an assessment for 2021 BCB410H:
Applied Bioinformatics, University of Toronto, Toronto, CANADA.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(Cov2Comparator)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
