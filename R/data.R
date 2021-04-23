#' T1D Application data
#'
#' A dataframe containing the chromosome (chrom), base pair position (pos),
#' rsid, major (other_allele) and minor (effect_allele) allele, T1D GWAS p-value (p),
#' binary chromatin accessibility measure (DGF_ENCODE), H3K27ac fold change value
#' in naive CD4+ T helper cells (Th_H3K27ac), RA GWAS p-value (RA_p), ldak weight, and minor allele frequency (maf) for 114,641 SNPs in the T1D GWAS (https://www.nature.com/articles/ng.3245)
#'
#' Minor allele frequency is taken from data by a study by Cooper et al.
#' (https://doi.org/10.1101/120022) where available. Missing values were replaced
#' by drawing samples from the empirical distribution of MAFs 
#'
#' @format A data frame with 114641 rows and 11 variables:
#' 
"T1D_df"
