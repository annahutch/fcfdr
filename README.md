
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fcfdr

<img src="man/figures/logo.png" align="right" />

`fcfdr` is an all-encompassing R package to implement the cFDR approach
for a variety of auxiliary covariates. As input it requires GWAS
p-values for a trait of interest and auxiliary data values. It then
outputs “v-values” which can be interpreted as GWAS p-values that have
been adjusted using the auxiliary data values.

If you have any questions please do not hesitate to contact me:
`annahutchinson1995@gmail.com`

Webpage: <https://annahutch.github.io/fcfdr/>

GitHub repository: <https://github.com/annahutch/fcfdr>

------------------------------------------------------------------------

## Manuscripts

**Flexible cFDR:**

Hutchinson A, Reales G, Willis T, Wallace C (2021) Leveraging auxiliary
data from arbitrary distributions to boost GWAS discovery with Flexible
cFDR. PLoS Genet 17(10): e1009853.
<https://doi.org/10.1371/journal.pgen.1009853>

**Binary cFDR:**

Hutchinson A, Liley J, Wallace C (2021) fcfdr: an R package to leverage
continuous and binary functional genomic data in GWAS. bioRxiv
2021.10.21.465274. <https://doi.org/10.1101/2021.10.21.465274>

------------------------------------------------------------------------

## Installation

You can install `fcfdr` from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("annahutch/fcfdr")
```

------------------------------------------------------------------------

## Examples

See the vignettes for examples of usage.

------------------------------------------------------------------------

## Brief history of cFDR

The cFDR framework was first introduced by [Andreassen and colleagues in
2013](https://doi.org/10.1371/journal.pgen.1003455). The group also
applied the approach to uncover new genetic associations for [blood
pressure](https://doi.org/10.1161/HYPERTENSIONAHA.113.02077) and
[multiple sclerosis](https://doi.org/10.1038/mp.2013.195) by leveraging
GWAS test statistics for related traits.

In 2015, [Liley and
Wallace](https://doi.org/10.1371/journal.pgen.1004926) showed that the
conventional cFDR approach does not control the frequentist FDR and
described an extension of the approach to ensure that FDR is controlled,
however their approach is rather conservative. In 2021, [Liley and
Wallace](https://doi.org/10.1002/bimj.201900254) described an elegant
cFDR framework that takes as input GWAS p-values for related traits and
returns “v-values” which can be interpreted as GWAS p-values that have
been adjusted using the auxiliary related trait data, and thus can be
used directly in any error-rate controlling procedure
(e.g. Benjamini-Hochberg).

Up until this point, the cFDR framework was designed for a very specific
setting, that is to increase GWAS discovery (in the “principal trait”)
by leveraging GWAS test statistics from a genetically related
(“conditional”) trait. In 2021 we described an extension of the cFDR
framework called Flexible cFDR that allows for auxiliary covariates
sampled from arbitrary continuous distributions to be leveraged with
GWAS test statistics - thus enabling broader applicability. We also
describe an extension for binary covariates, called Binary cFDR. In
practice, Flexible cFDR and Binary cFDR can be used to leverage
functional genomic data with GWAS results to boost power for GWAS
discovery whilst controlling the FDR.
