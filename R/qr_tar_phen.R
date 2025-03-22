#' A phenotype data file from the target dataset, containing the quantitative outcome treated as an exposure for the reverse direction.
#'
#' This dataset contains covariate information of individuals in the target dataset,
#' including various confounders.
#'
#' @format A data frame with 800 rows and 19 columns:
#' \describe{
#'   \item{Family ID}{Family identification number}
#'   \item{Individual ID}{Individual identification number}
#'   \item{quantitative Covariate}{A quantitative covariate}
#'   \item{Confounder 1}{Confounder 1 description}
#'   ...
#'   \item{Confounder 16}{Confounder 16 description}
#' }
#'
#' @usage data(qr_tar_phen)
#' @docType data
#' @keywords datasets
#' @name qr_tar_phen
NULL