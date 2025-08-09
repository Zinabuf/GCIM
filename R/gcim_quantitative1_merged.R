# R/gcim_quantitative.R

#' GCIM Analysis for Quantitative Outcomes
#'
#' Performs linear regression analysis for the Genetic Causality Inference Model (GCIM) 
#' with quantitative outcomes, without including the main effect of the exposure's PRS 
#' from own PRS data.
#'
#' @param qp_tar_phen File path to the phenotype and PRS data (tab-delimited). 
#' Must contain at least:
#'   FID, IID, Outcome, Add_PRS, Int_PRS, Cov_PRS, Covariate_Pheno.
#' Additional columns will be treated as confounders.
#' @param verbose Logical; whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \item{model}{The fitted `lm` model object.}
#' \item{model_summary}{Summary of the fitted model.}
#' \item{merged_data}{Data frame used for model fitting.}
#' \item{formula}{The model formula used.}
#' \item{sample_size}{Number of observations used in the analysis.}
#' \item{variables}{Names of the variables in the model.}
#' \item{r_squared}{R-squared value of the model.}
#' \item{adj_r_squared}{Adjusted R-squared value of the model.}
#' \item{residual_se}{Residual standard error of the model.}
#' \item{f_statistic}{F-statistic for the overall model.}
#'
#' @export
gcim_q1_merged <- function(qp_tar_phen, verbose = TRUE) {
  
  # Load phenotype & PRS data
  if (verbose) cat("Loading phenotype and PRS data...\n")
  merged_data <- utils::read.table(qp_tar_phen, header = TRUE)
  
  # Required columns
  required_cols <- c("FID", "IID", "Outcome", "Add_PRS", "Int_PRS", "Cov_PRS", "Covariate_Pheno")
  missing_cols <- setdiff(required_cols, colnames(merged_data))
  
  if (length(missing_cols) > 0) {
    stop("Error: Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Ensure outcome is numeric
  merged_data$Outcome <- as.numeric(merged_data$Outcome)
  
  # Check for valid outcome
  if (all(is.na(merged_data$Outcome))) {
    stop("Error: All outcome values are missing (NA)")
  }
  
  if (verbose) {
    cat(sprintf(
      "Outcome summary: Mean = %.3f, SD = %.3f, Range = [%.3f, %.3f]\n",
      mean(merged_data$Outcome, na.rm = TRUE),
      stats::sd(merged_data$Outcome, na.rm = TRUE),
      min(merged_data$Outcome, na.rm = TRUE),
      max(merged_data$Outcome, na.rm = TRUE)
    ))
  }
  
  # Build formula
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Cov_PRS + Int_PRS:Cov_PRS"
  
  # Add confounders if available
  confounder_vars <- setdiff(colnames(merged_data), required_cols)
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }
  
  model_formula <- stats::as.formula(formula_str)
  
  if (verbose) {
    cat("Fitting linear regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")
  }
  
  # Fit model
  result <- tryCatch({
    model <- stats::lm(model_formula, data = merged_data)
    model_summary <- summary(model)
    
    list(
      model = model,
      model_summary = model_summary,
      merged_data = merged_data,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      r_squared = model_summary$r.squared,
      adj_r_squared = model_summary$adj.r.squared,
      residual_se = model_summary$sigma,
      f_statistic = model_summary$fstatistic
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