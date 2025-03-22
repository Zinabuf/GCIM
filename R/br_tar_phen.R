#' A phenotype data file from target dataset, containing a binary outcome treated as a binary exposure in the reverse causal directions.
#'
#' This dataset contains phenotypic data for the covariate information of individuals in the target dataset,
#'  including various confounders.
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
#' @usage data(br_tar_phen)
#' @docType data
#' @keywords datasets
#' @name br_tar_phen
NULL