#' PLINK .map File for SNP Information
#'
#' This dataset follows the standard PLINK .map file format and contains genetic marker 
#' information, including chromosome, SNP ID, position in centimorgans, and base-pair 
#' coordinates.
#'
#' @format A data frame with 4 columns:
#' \describe{
#'   \item{Column 1}{Chromosome ID (integer)}
#'   \item{Column 2}{SNP ID (character)}
#'   \item{Column 3}{Position in centimorgans (numeric; typically 0 if unknown)}
#'   \item{Column 4}{Base-pair coordinate (integer; position on the chromosome in base pairs)}
#' }
#'
#' @usage data(dis_snp.map)
#' @docType data
#' @keywords datasets
#' @name dis_snp.map
NULL