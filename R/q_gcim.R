 # R/gcim_quantitative.R

#' Perform regression analysis for GCIM with continuous outcome
#'
#' This function performs linear regression analysis for Gene-Covariate Interaction 
#' Modeling (GCIM) with continuous outcomes using Polygenic Risk Scores (PRS).
#'
#' @param qp_tar_phen File path for the target phenotype data (FID, IID, Outcome format)
#' @param qp_tar_cov File path for the target covariate data (FID, IID, Covariate, Confounders format)
#' @param Add_PRS Either file path or data frame for additive PRS values
#' @param Int_PRS Either file path or data frame for interaction PRS values
#' @param Cov_PRS Either file path or data frame for covariate PRS values
#' @return List containing model summary and diagnostic information
#' @export
#' @examples
#' \dontrun{
#' result <- gcim_q("phenotype.txt", "covariates.txt", 
#'                          Add_PRS, Int_PRS, Cov_PRS)
#' }
gcim_q <- function(qp_tar_phen, qp_tar_cov, Add_PRS, Int_PRS, Cov_PRS){
  
  cat("Loading data files for quantitative outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- read.table(qp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is numeric
  outcome_data$Outcome <- as.numeric(outcome_data$Outcome)

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
  covariate_data <- read.table(qp_tar_cov, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  
  # Create column names for covariate data
  if(ncol(covariate_data) >= 3) {
    cov_colnames <- c("FID", "IID", "Covariate_Pheno")
    
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
  
    cat(sprintf("Final dataset contains %d observations with %d variables.\n", 
                nrow(merged_data), ncol(merged_data)))
  
  # Identify confounder variables (if any)
  confounder_vars <- colnames(merged_data)[grepl("^Conf_", colnames(merged_data))]
  
  # Build formula string
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS"
  
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }
  
  model_formula <- as.formula(formula_str)
  
    cat("Fitting linear regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")
  
  # Fit model
  tryCatch({
    model <- lm(model_formula, data = merged_data)
    model_summary <- summary(model)
    
    return(list(
      model_summary = model_summary,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data)
))
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    return(list(
      error = e$message,
      merged_data = merged_data,
      formula = model_formula
    ))
  })
}