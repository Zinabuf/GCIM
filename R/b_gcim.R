# R/gcim_binary.R

#' Perform regression analysis for GCIM with binary outcome
#'
#' This function performs logistic regression analysis for Gene-Covariate Interaction 
#' Modeling (GCIM) with binary outcomes using Polygenic Risk Scores (PRS).
#'
#' @param bp_tar_phen File path for the target phenotype data (FID, IID, Outcome format)
#' @param bp_tar_cov File path for the target covariate data (FID, IID, Covariate, Confounders format)
#' @param Add_PRS Either file path or data frame for additive PRS values
#' @param Int_PRS Either file path or data frame for interaction PRS values  
#' @param Cov_PRS Either file path or data frame for covariate PRS values
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @return List containing model summary and diagnostic information
#' @export
#' @examples
#' \dontrun{
#' # Using file paths
#' result <- gcim_b("phenotype.txt", "covariates.txt", 
#'                      "add_prs.txt", "int_prs.txt", "cov_prs.txt")
#' 
#' # Using data frames (from GxEprs output)
#' result <- gcim_b("phenotype.txt", "covariates.txt", 
#'                      add_prs_df, int_prs_df, cov_prs_df)
#' }
gcim_b <- function(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, verbose = TRUE) {
  
  if(verbose) cat("Loading data files for binary outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- read.table(bp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is binary (0/1)
  if(!all(outcome_data$Outcome %in% c(0, 1, NA))) {
    if(verbose) cat("Warning: Outcome variable contains non-binary values. Converting to binary.\n")
    outcome_data$Outcome <- as.numeric(as.factor(outcome_data$Outcome)) - 1
  }
  
  # Handle PRS data (could be file paths or data frames)
  # If q, r, p are file paths:
add_prs <- q
int_prs <- r
cov_prs <- p

# Rename 3rd column to standard expected names
colnames(add_prs)[1:3] <- c("FID", "IID", "Add_PRS")
colnames(int_prs)[1:3] <- c("FID", "IID", "Int_PRS")
colnames(cov_prs)[1:3] <- c("FID", "IID", "Cov_PRS")
  
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
  
  # Check if there are both cases and controls
  outcome_table <- table(merged_data$Outcome)
  if(length(outcome_table) < 2) {
    stop("Error: Binary outcome must have both cases (1) and controls (0)")
  }
  
  if(verbose) {
    cat(sprintf("Final dataset contains %d observations (%d cases, %d controls) with %d variables.\n", 
                nrow(merged_data), 
                sum(merged_data$Outcome == 1, na.rm = TRUE),
                sum(merged_data$Outcome == 0, na.rm = TRUE),
                ncol(merged_data)))
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
    cat("Fitting logistic regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")
  }

  # Fit logistic regression model
  tryCatch({
    model <- glm(model_formula, data = merged_data, family = binomial(link = "logit"))
    model_summary <- summary(model)
    
    # Calculate pseudo R-squared (McFadden's R-squared)
    null_deviance <- model$null.deviance
    residual_deviance <- model$deviance
    pseudo_r_squared <- 1 - (residual_deviance / null_deviance)
    
    return(list(
      model_summary = model_summary,
      merged_data = merged_data,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      pseudo_r_squared = pseudo_r_squared,
      null_deviance = null_deviance,
      residual_deviance = residual_deviance,
      aic = model$aic
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