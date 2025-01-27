#' A phenotype data from the discovery dataset containing binary phenotype trated as exposure for the proposed directions of causations. 
#'
#' This dataset contains covariate information for individuals in the discovery dataset,
#' including various confounders for a binary outcome.
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
#' @usage data(bp_dis_cov)
#' @docType data
#' @keywords datasets
#' @name bp_dis_cov
NULL