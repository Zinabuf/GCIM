# R/gcim_binary.R

#' Perform regression analysis for GCIM with binary outcome
#'
#' This function performs logistic regression analysis for Genetic causality inference model (GCIM) with binary outcomes using Polygenic Risk Scores (PRS).
#'
#' @param bp_tar_phen File path for the target phenotype data (FID, IID, Outcome format)
#' @param bp_tar_cov File path for the target covariate data (FID, IID, Covariate, Confounders format)
#' @param Add_PRS  data frame for additive PRS values
#' @param Int_PRS  data frame for interaction PRS values  
#' @param Cov_PRS  data frame for covariate PRS values
#' @return List containing model summary and diagnostic information
#' @export
#' @examples
#' \dontrun{
#' # Using data frames (from GxEprs output)
#' result <- gcim_b("phenotype.txt", "covariates.txt", 
#'                      Add_PRS, Int_PRS, Cov_PRS)
#' }
gcim_b <- function(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS) {
  
   cat("Loading data files for binary outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- read.table(bp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is binary (0/1)
  if(!all(outcome_data$Outcome %in% c(0, 1, NA))) {
     cat("Warning: Outcome variable contains non-binary values. Converting to binary.\n")
    outcome_data$Outcome <- as.numeric(as.factor(outcome_data$Outcome)) - 1
  }
  
  # Create Add_PRS from your 'q' data (additive PRS)
add_prs <- data.frame(
  FID = q$FID,
  IID = q$IID,
  Add_PRS = q[,3])  # Replace [,3] with actual column name if needed

# Create Int_PRS from your 'r' data (interaction PRS)  
int_prs <- data.frame(
  FID = r$FID,
  IID = r$IID,
  Int_PRS = r[,3])  # Replace [,3] with actual column name if needed

# Create Cov_PRS from your 'p' data (covariate PRS)
cov_prs <- data.frame(
  FID = p$FID,
  IID = p$IID,
  Cov_PRS = p[,3])  # Replace [,3] with actual column name if needed

# Scale the PRS scores
add_prs[,3] <- scale(add_prs[,3])
int_prs[,3] <- scale(int_prs[,3])
cov_prs[,3] <- scale(cov_prs[,3])
  
  # Load covariate data
  covariate_data <- read.table(bp_tar_cov, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  
  # Create column names for covariate data
  if(ncol(covariate_data) >= 3) {
    cov_colnames <- c("FID", "IID", "Covariate_Pheno")
    
    # Add confounder column names for remaining columns
    if(ncol(covariate_data) > 3) {
      n_confounders <- ncol(covariate_data) - 3
      conf_names <- paste0("Conf_", 1:n_confounders)
      cov_colnames <- c(cov_colnames, conf_names)
    }
    
    colnames(covariate_data) <- cov_colnames
  }
  
  # Merge data
   cat("Merging data files by FID...\n")
  
# Start with the outcome data
merged_data <- outcome_data

# Merge everything using both FID and IID
merged_data <- merge(merged_data, add_prs[, c("FID", "IID", "Add_PRS")], by = c("FID", "IID"), all.x = TRUE)
merged_data <- merge(merged_data, int_prs[, c("FID", "IID", "Int_PRS")], by = c("FID", "IID"), all.x = TRUE)
merged_data <- merge(merged_data, cov_prs[, c("FID", "IID", "Cov_PRS")], by = c("FID", "IID"), all.x = TRUE)
merged_data <- merge(merged_data, covariate_data, by = c("FID", "IID"), all.x = TRUE) 
      
  # Check if there are both cases and controls
  outcome_table <- table(merged_data$Outcome)
  if(length(outcome_table) < 2) {
    stop("Error: Binary outcome must have both cases (1) and controls (0)")
  }
  
    cat(sprintf("Final dataset contains %d observations (%d cases, %d controls) with %d variables.\n", 
                nrow(merged_data), 
                sum(merged_data$Outcome == 1, na.rm = TRUE),
                sum(merged_data$Outcome == 0, na.rm = TRUE),
                ncol(merged_data)))
  
  # Identify confounder variables (if any)
  confounder_vars <- colnames(merged_data)[grepl("^Conf_", colnames(merged_data))]
  
  # Build formula string
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS"
  
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }

  model_formula <- as.formula(formula_str)

    cat("Fitting logistic regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")

 # Fit logistic regression model
tryCatch({
  model <- glm(model_formula, data = merged_data, family = binomial(link = "logit"))
  model_summary <- summary(model)

  null_deviance <- model$null.deviance
  residual_deviance <- model$deviance

  return(list(
    model_summary = model_summary,
    merged_data = merged_data,
    formula = model_formula,
    sample_size = nrow(merged_data),
    variables = colnames(merged_data),
    null_deviance = null_deviance,
    residual_deviance = residual_deviance,
    aic = model$aic
  ))
}, error = function(e) {
  warning(paste("Model fitting failed:", e$message))
  return(list(
    error = e$message,
    merged_data = merged_data,
    formula = model_formula
  ))
})
}