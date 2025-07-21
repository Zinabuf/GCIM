 # R/gcim_continuous.R

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
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @return List containing model summary and diagnostic information
#' @export
#' @examples
#' \dontrun{
#' result <- gcim_q("phenotype.txt", "covariates.txt", 
#'                          add_prs, int_prs, cov_prs)
#' }
gcim_q <- function(qp_tar_phen, qp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, verbose = TRUE){
  
  if(verbose) cat("Loading data files for continuous outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- read.table(qp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is numeric
  outcome_data$Outcome <- as.numeric(outcome_data$Outcome)
  
  # Handle PRS data (could be file paths or data frames)
# If q, r, p are file paths:
add_prs <- q
int_prs <- r
cov_prs <- p

# Rename 3rd column to standard expected names
colnames(add_prs)[3] <- "Add_PRS"
colnames(int_prs)[3] <- "Int_PRS"
colnames(cov_prs)[3] <- "Cov_PRS"
  
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
  if(verbose) cat("Merging data files by FID...\n")
  
  merged_data <- outcome_data
  merged_data <- merge(merged_data, add_prs[, c("FID", "Add_PRS")], by = "FID", all.x = TRUE)
  merged_data <- merge(merged_data, int_prs[, c("FID", "Int_PRS")], by = "FID", all.x = TRUE)
  merged_data <- merge(merged_data, cov_prs[, c("FID", "Cov_PRS")], by = "FID", all.x = TRUE)
  merged_data <- merge(merged_data, covariate_data, by = "FID", all.x = TRUE)
  
  # Remove rows with missing essential variables
  essential_vars <- c("Outcome", "Add_PRS", "Int_PRS", "Cov_PRS", "Covariate_Pheno")
  complete_cases <- complete.cases(merged_data[, essential_vars])
  merged_data <- merged_data[complete_cases, ]
  
  if(verbose) {
    cat(sprintf("Final dataset contains %d observations with %d variables.\n", 
                nrow(merged_data), ncol(merged_data)))
  }
  
  # Identify confounder variables (if any)
  confounder_vars <- colnames(merged_data)[grepl("^Conf_", colnames(merged_data))]
  
  # Build formula string
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS"
  
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }
  
  model_formula <- as.formula(formula_str)
  
  if(verbose) {
    cat("Fitting linear regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")
  }
  
  # Fit model
  tryCatch({
    model <- lm(model_formula, data = merged_data)
    model_summary <- summary(model)
    
    return(list(
      model_summary = model_summary,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      r_squared = model_summary$r.squared,
      adj_r_squared = model_summary$adj.r.squared
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