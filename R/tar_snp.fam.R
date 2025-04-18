#' PLINK .fam file for family and phenotype information
#'
#' This dataset follows the standard PLINK .fam format and contains family structure 
#' and phenotype information for individuals in the target dataset.
#'
#' @format A data frame with 6 columns:
#' \describe{
#'   \item{Column 1}{Family ID (character or integer)}
#'   \item{Column 2}{Individual ID (character or integer)}
#'   \item{Column 3}{Father's ID (character or integer; 0 if missing)}
#'   \item{Column 4}{Mother's ID (character or integer; 0 if missing)}
#'   \item{Column 5}{Sex (integer; 1 = male, 2 = female, 0 = unknown)}
#'   \item{Column 6}{Phenotype value (numeric; 1 = control, 2 = case, -9 = missing)}
#' }
#'
#' @usage data(tar_snp.fam)
#' @docType data
#' @keywords datasets
#' @name tar_snp.fam
NULL