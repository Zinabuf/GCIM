#' Perform regression analysis for GCIM with quantitative outcome.
#'
#' @param qp_tar_phen File path for the target phenotype data.
#' @param qp_tar_cov File path for the target covariate data.
#' @param Add_PRS File path for scaled additive PRS values.
#' @param Int_PRS File path for scaled interaction PRS values.
#' @param Cov_PRS File path for scaled covariate PRS values.
#' @param confounders Optional data frame of additional confounders (not used if covariates file contains confounders).
#' @return Summary of the regression model.
#' @export
gcim_q <- function(qp_tar_phen, qp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders = NULL) {
  temp_dir <- tempdir()  
  
  # Load all data files
  cat("Loading data files...\n")
  
  # Load phenotype data
  outcome_data <- read.table(qp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Load PRS data files
  prs_add <- read.table(Add_PRS, header = TRUE, stringsAsFactors = FALSE)
    colnames(prs_add)[1:3] <- c("FID", "IID", "Add_PRS")
  
  prs_int <- read.table(Int_PRS, header = TRUE, stringsAsFactors = FALSE)
    colnames(prs_int)[1:3] <- c("FID", "IID", "Int_PRS")
  
  prs_cov <- read.table(Cov_PRS, header = TRUE, stringsAsFactors = FALSE)
    colnames(prs_cov)[1:3] <- c("FID", "IID", "Cov_PRS")
  
  # Load covariate data
  covariate_data <- read.table(qp_tar_cov, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  
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
  
  # Start merging by FID
  cat("Merging data files by FID...\n")
  
  # Start with outcome data
  merged_data <- outcome_data
  
  # Merge PRS data
  merged_data <- merge(merged_data, prs_add[, c("FID", "Add_PRS")], by = "FID", all.x = TRUE)
  merged_data <- merge(merged_data, prs_int[, c("FID", "Int_PRS")], by = "FID", all.x = TRUE)
  merged_data <- merge(merged_data, prs_cov[, c("FID", "Cov_PRS")], by = "FID", all.x = TRUE)
  
  # Merge covariate data
  merged_data <- merge(merged_data, covariate_data, by = "FID", all.x = TRUE)
  
  # Add external confounders if provided
  if(!is.null(confounders)) {
    if("FID" %in% colnames(confounders)) {
      merged_data <- merge(merged_data, confounders, by = "FID", all.x = TRUE)
    } else {
      cat("Warning: External confounders do not contain FID column. Skipping external confounders.\n")
    }
  }
  
  # Remove rows with missing essential variables
  essential_vars <- c("Outcome", "Add_PRS", "Int_PRS", "Cov_PRS", "Covariate_Pheno")
  complete_cases <- complete.cases(merged_data[, essential_vars])
  merged_data <- merged_data[complete_cases, ]
  
  cat(sprintf("Final dataset contains %d observations with %d variables.\n", 
              nrow(merged_data), ncol(merged_data)))
  
  # Identify confounder columns (all columns except FID, IID, and essential variables)
  exclude_cols <- c("FID", "IID", "Outcome", "Add_PRS", "Int_PRS", "Cov_PRS", "Covariate_Pheno")
  confounder_cols <- setdiff(colnames(merged_data), exclude_cols)
  
  # Build regression formula
  if(length(confounder_cols) > 0) {
    # Include confounders in the model
    confounders_formula <- paste(confounder_cols, collapse = " + ")
    model_formula <- as.formula(paste("Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS +", 
                                      confounders_formula))
  } else {
    # No confounders available
    model_formula <- as.formula("Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS")
  }
  
  cat("Fitting regression model...\n")
  cat("Model formula:", deparse(model_formula), "\n")
  
  # Fit the regression model
  tryCatch({
    model <- lm(model_formula, data = merged_data)
    
    # Print data summary
    cat("\nData Summary:\n")
    cat("Variables in final dataset:\n")
    print(colnames(merged_data))
    cat("\nSample size:", nrow(merged_data), "\n")
    
    # Return both model summary and merged data
    return(list(
      model_summary = summary(model),
      merged_data = merged_data,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data)
    ))
    
  }, error = function(e) {
    cat("Error fitting model:", e$message, "\n")
    return(list(
      error = e$message,
      merged_data = merged_data,
      formula = model_formula
    ))
  })
}
