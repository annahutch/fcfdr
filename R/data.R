#' Data for T1D application
#'
#' A data.frame containing the rsID, chromosome (CHR19) and base pair position (BP19) in hg19,
#' reference allele (REF), alternative allele (ALLT), type 1 diabetes GWAS p-value (T1D_pval),
#' minor allele frequency (MAF), LDAK weight (LDAK_weight), rheumatoid arthritis GWAS p-value (RA_pval),
#' binary regulatory factor binding site overlap (DGF), average H3K27ac fold change value
#' in T1D-relevant cell types (H3K27ac) for 113,543 SNPs in the T1D GWAS (https://www.nature.com/articles/ng.3245)
#'
#' Minor allele frequencies estimated from the CEU sub-population samples 
#' in the 1000 Genomes Project Phase 3 data set. Missing values were replaced
#' by drawing samples from the empirical distribution of MAFs 
#'
#' @format A data frame with 113543 rows and 11 variables:
#' 
"T1D_application_data"




# Suppress R CMD check notes
#' @importFrom graphics hist
#' @importFrom graphics image
#' @importFrom graphics lines
#' @importFrom dplyr transmute
#' @importFrom data.table as.xts.data.table
NULL