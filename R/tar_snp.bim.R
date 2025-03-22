#' PLINK .bim file for SNP information
#'
#' This dataset follows the standard PLINK .bim format and contains SNP-related information
#' for the target dataset, including chromosome ID, SNP ID, and allele details.
#'
#' @format A data frame with 6 columns:
#' \describe{
#'   \item{Column 1}{Chromosome ID (integer)}
#'   \item{Column 2}{SNP ID (character)}
#'   \item{Column 3}{Position in centimorgans (numeric)}
#'   \item{Column 4}{Base-pair coordinate (integer)}
#'   \item{Column 5}{Minor Allele (character)}
#'   \item{Column 6}{Reference Allele (character)}
#' }
#'
#' @usage data(tar_snp.bim)
#' @docType data
#' @keywords datasets
#' @name tar_snp.bim
NULL