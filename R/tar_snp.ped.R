#' PLINK .ped file for genotype and sample information
#'
#' This dataset follows the standard PLINK .ped file format. It contains family, 
#' individual, parental, sex, phenotype, and genotype information for each sample.
#'
#' @format A data frame with 6 required columns, followed by genotype data:
#' \describe{
#'   \item{Column 1}{Family ID (character)}
#'   \item{Column 2}{Individual ID (character)}
#'   \item{Column 3}{Paternal ID (0 if unknown, character)}
#'   \item{Column 4}{Maternal ID (0 if unknown, character)}
#'   \item{Column 5}{Sex (1 = male, 2 = female, 0 = unknown)}
#'   \item{Column 6}{Phenotype value (1 = control, 2 = case, -9 = missing)}
#'   \item{Columns 7+}{Genotype data in pairs of alleles for each SNP 
#'       (e.g., A T, C, G, coded as separate columns).}
#' }
#'
#' @usage data(tar_snp.ped)
#' @docType data
#' @keywords datasets
#' @name tar_snp.ped
NULL