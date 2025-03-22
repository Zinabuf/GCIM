#'  A phenotype data from the discovery dataset containing quantitative phenotype trated as exposure for the proposed causal directions. 
#'
#' This dataset contains covariate information for individuals in the discovery dataset,
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
#' @usage data(qp_dis_cov)
#' @docType data
#' @keywords datasets
#' @name qp_dis_cov
NULL