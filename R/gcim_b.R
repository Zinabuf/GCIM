#' Perform regression analysis for GCIM with binary outcome.
#'
#' @param bp_tar_phen File path for the target phenotype data.
#' @param bp_tar_cov File path for the target covariate data.
#' @param Additive Scaled additive PRS values.
#' @param Interaction Scaled interaction PRS values.
#' @param Covariate Scaled covariate PRS values.
#' @param Confounders Data frame of additional confounders.
#' @return Summary of the regression model.
#' @export
gcim_b <- function(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders) {
  temp_dir <- tempdir()  
  # Load phenotype and covariate data
  outcome_bp_data <- read.table(bp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  prs_add <- read.table(file.path(temp_dir, "prs_add_scaled.txt"), header = TRUE, stringsAsFactors = FALSE)
  prs_int <- read.table(file.path(temp_dir, "prs_int_scaled.txt"), header = TRUE, stringsAsFactors = FALSE)
  prs_cov <- read.table(file.path(temp_dir, "prs_cov_scaled.txt"), header = TRUE, stringsAsFactors = FALSE)
  covariate_bp_data <- read.table(bp_tar_cov, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)

  # Prepare regression data
  regression_data <- data.frame(
    Outcome = outcome_bp_data[, 3],  # Assuming 3rd column is the outcome
    Add_PRS = prs_add[, 3],   # Adjust column index if necessary
    Int_PRS = prs_int[, 3],
    Cov_PRS = prs_cov[, 3],
    Covariate_Pheno = covariate_bp_data[, 3] # Adjust column index as needed
  )

  # Add confounders dynamically based on the input `conf_`
  confounders <- covariate_bp_data[, 4:19]
  colnames(confounders) <- paste0("Conf_", seq_along(confounders))  # Rename confounders

  # Combine confounders with regression data
  regression_data <- cbind(regression_data, confounders)

  # Fit the regression model using all variables
  model_formula <- as.formula(paste("Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS +", paste(names(confounders), collapse = " + ")))
  model <- glm(model_formula, family = binomial(), data = regression_data)

  # Return model summary
  return(summary(model))
}