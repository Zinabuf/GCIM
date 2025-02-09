
#' A phenotype data file from the target dataset, containing the binary exposure treated as an oucome for the reverse direction of causation 
#'
#' This dataset contains phenotype data for the outcome information for individuals in the target dataset,
#'  for a binary outcome as a binary outcome in the reverse directions of causations. 
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
#' @usage data(br_tar_cov)
#' @docType data
#' @keywords datasets
#' @name br_tar_cov
NULL