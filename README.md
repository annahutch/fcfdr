
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcfdr

<img src="man/figures/logo.png" align="right" />

`fcfdr` is an R package that implements cFDR to generate GWAS
\(p\)-values that have been adjusted using auxiliary covariates. It
supports auxiliary covariates from any arbitrary distribution, enabling
researchers to leverage a wide variety of auxiliary data types with GWAS
\(p\)-values to boost power for discovery.

If you have any questions please do not hesitate to contact me:
`anna.hutchinson@mrc-bsu.cam.ac.uk`

Webpage: <https://annahutch.github.io/fcfdr/>

-----

## Installation

You can install `fcfdr` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("annahutch/fcfdr")
```

-----

## Examples

See the vignettes for examples of usage.

-----

## Abstract

Genome-wide association studies (GWAS) have identified thousands of
genetic variants that are associated with complex traits. However, a
stringent significance threshold is required to identify robust genetic
associations. Leveraging relevant auxiliary covariates has the potential
to boost statistical power to exceed the significance threshold.
Particularly, abundant pleiotropy and the non-random distribution of
SNPs across various functional categories suggests that leveraging GWAS
test statistics from related traits and/or functional genomic data may
boost GWAS discovery. While type 1 error rate control has become
standard in GWAS, control of the false discovery rate can be a more
powerful approach. The conditional false discovery rate (cFDR) extends
the standard FDR framework by conditioning on auxiliary data to call
significant associations, but current implementations are restricted to
auxiliary data satisfying specific parametric distributions, typically
GWAS \(p\)-values for related traits. We relax these distributional
assumptions, enabling an extension of the cFDR framework that supports
auxiliary covariates from arbitrary continuous distributions (“Flexible
cFDR”). Our method can be applied iteratively, thereby supporting
multi-dimensional covariate data. Through simulations we show that
Flexible cFDR increases sensitivity whilst controlling FDR after one or
several iterations. We further demonstrate its practical potential
through application to an asthma GWAS, leveraging various functional
genomic data to find additional genetic associations for asthma, which
we validate in the larger, independent, UK Biobank data resource.
