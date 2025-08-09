# R/gcim_b1.R
#' Perform regression analysis for genetic causality inference model(GCIM) with binary outcome
#'
#' This function performs logistic regression analysis for GCIM 
#' with binary outcomes using the main effects of Polygenic Risk Scores (PRS) of the exposure.
#' It can read PRS files from temporary directories or use provided PRS data objects.
#' from own PRS data.
#'
#' @param bp_tar_phen File path to the phenotype and PRS data (tab-delimited). 
#' Must contain at least:
#'   FID, IID, Outcome, Add_PRS, Int_PRS, Cov_PRS, Covariate_Pheno.
#' Additional columns will be treated as confounders.
#' @param verbose Logical; whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \item{model}{The fitted `glm` model object.}
#' \item{model_summary}{Summary of the fitted model.}
#' \item{merged_data}{Data frame used for model fitting.}
#' \item{formula}{The model formula used.}
#' \item{sample_size}{Number of observations used in the analysis.}
#' \item{variables}{Names of the variables in the model.}
#' \item{pseudo_r_squared}{Pseudo R-squared value of the model.}
#' \item{null_deviance}{Null deviance of the model.}
#' \item{residual_deviance}{Residual deviance of the model.}
#' \item{aic}{AIC of the model.}
#'
#' @export
gcim_b1_merged <- function(bp_tar_phen, verbose = TRUE) {
  
  # Load phenotype & PRS data
  if (verbose) cat("Loading phenotype and PRS data...\n")
  merged_data <- utils::read.table(bp_tar_phen, header = TRUE)
  
  # Required columns
  required_cols <- c("FID", "IID", "Outcome", "Add_PRS", "Int_PRS", "Cov_PRS", "Covariate_Pheno")
  missing_cols <- setdiff(required_cols, colnames(merged_data))
  
  # FIXED: Use merged_data instead of undefined outcome_data
  if(!all(merged_data$Outcome %in% c(0, 1, NA))) {
    cat("Warning: Outcome variable contains non-binary values. Converting to binary.\n")
  }
  
  if (length(missing_cols) > 0) {
    stop("Error: Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check for valid outcome
  if (all(is.na(merged_data$Outcome))) {
    stop("Error: All outcome values are missing (NA)")
  }
  
  # Build formula (includes Cov_PRS main effect)
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Cov_PRS + Int_PRS:Cov_PRS"
  
  # Add confounders if available
  confounder_vars <- setdiff(colnames(merged_data), required_cols)
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }
  
  model_formula <- stats::as.formula(formula_str)
  
  if (verbose) {
    cat("Fitting logistic regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")
  }
  
  # Fit model
  result <- tryCatch({
    model <- stats::glm(model_formula, data = merged_data, family = stats::binomial(link = "logit"))
    model_summary <- summary(model)
    
    # FIXED: Define variables before using them in the list
    null_deviance <- model$null.deviance
    residual_deviance <- model$deviance
    pseudo_r_squared <- 1 - (residual_deviance / null_deviance)
    
    list(
      model = model,
      model_summary = model_summary,
      merged_data = merged_data,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      pseudo_r_squared = pseudo_r_squared,
      null_deviance = null_deviance,
      residual_deviance = residual_deviance,
      aic = model$aic
    )
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    list(
      error = e$message,
      merged_data = merged_data,
      formula = model_formula
    )
  })
  
  return(result)
}