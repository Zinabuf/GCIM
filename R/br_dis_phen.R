
#' A phenotype data file from the discovery dataset, containing the binary oucome treated as binary exposure for the reverse direction.
#'
#' This dataset contains phenotype information for individuals in the discovery dataset,
#' for a binary outcome treated as a binary exposure in the reverse direction.
#'
#' @format A data frame with 800 rows and 19 columns:
#' \describe{
#'   \item{Family ID}{Family identification number}
#'   \item{Individual ID}{Individual identification number}
#'   \item{Binary Covariate}{A binary covariate}
#'   \item{Confounder 1}{Confounder 1 description}
#'   ...
#'   \item{Confounder 16}{Confounder 16 description}
#' }
#'
#' @usage data(br_dis_phen)
#' @docType data
#' @keywords datasets
#' @name br_dis_phen
NULL