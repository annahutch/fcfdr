
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcfdr

<img src="man/figures/logo.png" align="right" />

`fcfdr` is an R package that implements cFDR to generate GWAS
\(p\)-values that have been adjusted using auxiliary covariates. It
supports auxiliary covariates from any arbitrary distribution, enabling
researchers to leverage a wide variety of auxiliary data types with GWAS
\(p\)-values.

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
genetic variants that are associated with complex traits. However,
erroneous associations become more likely when conducting many tests in
parallel, and consequently a stringent significance threshold is
required to identify robust genetic associations. Leveraging relevant
auxiliary data has the potential to boost the statistical power required
to exceed the significance threshold. Particularly, abundant pleiotropy
and the non-random distribution of GWAS SNPs across various functional
categories suggests that leveraging test statistics from related traits
and/or functional genomic data may boost GWAS discovery. The conditional
false discovery rate (cFDR) extends the standard FDR framework by
conditioning on auxiliary data to call significant associations. Current
methods estimate cFDR values using empirical cumulative distributions
and are restricted to leveraging auxiliary data such as \(p\)-values
measuring associations for related traits, which are expected to have a
known statistical distribution. Utilising bivariate kernel density
estimates (KDEs) rather than empirical estimates, we relax the
distributional assumptions, enabling an extension of the cFDR framework
that supports auxiliary data from any distribution (“Flexible cFDR”).
