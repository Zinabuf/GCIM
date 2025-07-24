# R/gcim_quantitative.R

#' Perform regression analysis for genetic causality inference model(GCIM) with quantitative outcome
#'
#' This function performs linear regression analysis for GCIM 
#' with quantitative outcomes using Polygenic Risk Scores (PRS).
#' It reads PRS files that were previously saved by PRS_quantitative function.
#'
#' @param qp_tar_phen File path for the target phenotype data (FID, IID, Outcome format)
#' @param qp_tar_cov File path for the target covariate data (FID, IID, Covariate, Confounders format)
#' @param Add_PRS Either file path or data frame for additive PRS values (fallback if files not found)
#' @param Int_PRS Either file path or data frame for interaction PRS values (fallback if files not found)
#' @param Cov_PRS Either file path or data frame for covariate PRS values (fallback if files not found)
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param scale_prs Logical, whether to scale PRS values (default: TRUE)
#' @param save_temp_files Logical, whether to save temporary PRS files (default: TRUE)
#' @return List containing model summary and diagnostic information
#' @export
#' @examples
#' \dontrun{
#' # After running PRS_quantitative functions to generate PRS files:
#' # add_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
#' # int_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe)  
#' # cov_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)
#' 
#' result <- gcim_q("Qphe_target.txt", "Qexp_target.txt", 
#'                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = cov_prs)
#' }
gcim_q <- function(qp_tar_phen, qp_tar_cov, Add_PRS, Int_PRS, Cov_PRS,
                   verbose = TRUE, scale_prs = TRUE, save_temp_files = TRUE) {
  
  if(verbose) cat("Loading data files for quantitative outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- utils::read.table(qp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is numeric
  outcome_data$Outcome <- as.numeric(outcome_data$Outcome)
  
  # Check for valid quantitative outcome
  if(all(is.na(outcome_data$Outcome))) {
    stop("Error: All outcome values are missing (NA)")
  }
  
  if(verbose) {
    cat(sprintf("Outcome summary: Mean = %.3f, SD = %.3f, Range = [%.3f, %.3f]\n",
                mean(outcome_data$Outcome, na.rm = TRUE),
                stats::sd(outcome_data$Outcome, na.rm = TRUE),
                min(outcome_data$Outcome, na.rm = TRUE),
                max(outcome_data$Outcome, na.rm = TRUE)))
  }
  
  # Read PRS files that were saved by PRS_quantitative function in temporary directories
  if(verbose) cat("Reading saved PRS files from temporary directories...\n")
  
  # Get the temporary directory path
  temp_dir <- tempdir()
  if(verbose) cat("Using temporary directory:", temp_dir, "\n")
  
  # Function to read PRS file with flexible naming
  read_prs_file <- function(prs_type, fallback_data = NULL) {
    # Try different possible file names
    possible_files <- c(
      file.path(temp_dir, paste0(prs_type, "_prs.sscore")),
      file.path(temp_dir, paste0(prs_type, ".sscore")), 
      file.path(temp_dir, "prs.sscore"),
      file.path(temp_dir, paste0(prs_type, "_prs.txt")),
      file.path(temp_dir, paste0(prs_type, ".txt"))
    )
    
    for(file_path in possible_files) {
      if(file.exists(file_path)) {
        if(verbose) cat("Found", prs_type, "PRS file:", basename(file_path), "\n")
        prs_data <- utils::read.table(file_path, header = FALSE)
        
        # Handle different file formats
        if(ncol(prs_data) >= 5) {
          # Standard plink .sscore format (FID, IID, PHENO1, CNT, SCORE1, SCORE2_AVG)
          prs_data <- prs_data[, c(1, 2, 5)]  # Columns 1, 2, 5 as per PRS_quantitative
        } else if(ncol(prs_data) >= 3) {
          # Simplified format (FID, IID, PRS)
          prs_data <- prs_data[, c(1, 2, 3)]
        } else {
          next  # Try next file
        }
        
        colnames(prs_data) <- c("FID", "IID", "PRS")
        return(prs_data)
      }
    }
    
    # If no file found, use fallback data
    if(!is.null(fallback_data)) {
      if(verbose) cat("No", prs_type, "PRS file found, using provided data object\n")
      # Ensure fallback data has correct format
      if(is.data.frame(fallback_data) && ncol(fallback_data) >= 3) {
        fallback_data <- fallback_data[, c(1, 2, 3)]
        colnames(fallback_data) <- c("FID", "IID", "PRS")
        return(fallback_data)
      } else if(is.character(fallback_data) && length(fallback_data) == 1 && file.exists(fallback_data)) {
        # If fallback_data is a file path
        prs_data <- utils::read.table(fallback_data, header = FALSE)
        if(ncol(prs_data) >= 3) {
          prs_data <- prs_data[, c(1, 2, 3)]
          colnames(prs_data) <- c("FID", "IID", "PRS")
          return(prs_data)
        }
      }
    }
    
    stop("No ", prs_type, " PRS file found and no valid fallback data provided")
  }
  
  # Read each PRS type
  if(verbose) cat("Reading additive PRS...\n")
  add_prs_raw <- read_prs_file("add", Add_PRS)
  
  if(verbose) cat("Reading interaction PRS...\n") 
  int_prs_raw <- read_prs_file("int", Int_PRS)
  
  if(verbose) cat("Reading covariate PRS...\n")
  cov_prs_raw <- read_prs_file("cov", Cov_PRS)
  
  # Save copies to standardized temporary files if requested
  temp_files <- NULL
  if(save_temp_files) {
    if(verbose) cat("Saving PRS copies to standardized temporary files...\n")
    
    # Define standardized temporary file paths  
    add_prs_file <- file.path(temp_dir, "add_prs.txt")
    int_prs_file <- file.path(temp_dir, "int_prs.txt")
    cov_prs_file <- file.path(temp_dir, "cov_prs.txt")
    
    # Save copies with headers for easier reading
    utils::write.table(add_prs_raw, add_prs_file, quote = FALSE, row.names = FALSE, sep = "\t")
    utils::write.table(int_prs_raw, int_prs_file, quote = FALSE, row.names = FALSE, sep = "\t")
    utils::write.table(cov_prs_raw, cov_prs_file, quote = FALSE, row.names = FALSE, sep = "\t")
    
    # Store file paths for return object
    temp_files <- c(add_prs_file, int_prs_file, cov_prs_file)
  }
  
  # Process the PRS data for analysis
  if(verbose) cat("Processing PRS data for analysis...\n")
  add_prs <- add_prs_raw
  int_prs <- int_prs_raw  
  cov_prs <- cov_prs_raw
  
  # Standardize column names for PRS data
  colnames(add_prs) <- c("FID", "IID", "Add_PRS")
  colnames(int_prs) <- c("FID", "IID", "Int_PRS")
  colnames(cov_prs) <- c("FID", "IID", "Cov_PRS")
  
  # Scale the PRS scores if requested
  if(scale_prs) {
    if(verbose) cat("Scaling PRS values (mean=0, sd=1)...\n")
    add_prs$Add_PRS <- as.numeric(scale(add_prs$Add_PRS)[,1])
    int_prs$Int_PRS <- as.numeric(scale(int_prs$Int_PRS)[,1])
    cov_prs$Cov_PRS <- as.numeric(scale(cov_prs$Cov_PRS)[,1])
  }
  
  if(verbose) {
    cat("PRS data summary:\n")
    cat("- Additive PRS: ", nrow(add_prs), " individuals\n")
    cat("- Interaction PRS: ", nrow(int_prs), " individuals\n") 
    cat("- Covariate PRS: ", nrow(cov_prs), " individuals\n")
  }
  
  # Load covariate data
  covariate_data <- utils::read.table(qp_tar_cov, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  
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
  if(verbose) cat("Merging data files by FID and IID...\n")
  
  # Start with the outcome data
  merged_data <- outcome_data
  
  # Merge everything using both FID and IID
  merged_data <- merge(merged_data, add_prs[, c("FID", "IID", "Add_PRS")], by = c("FID", "IID"), all.x = TRUE)
  merged_data <- merge(merged_data, int_prs[, c("FID", "IID", "Int_PRS")], by = c("FID", "IID"), all.x = TRUE)
  merged_data <- merge(merged_data, cov_prs[, c("FID", "IID", "Cov_PRS")], by = c("FID", "IID"), all.x = TRUE)
  merged_data <- merge(merged_data, covariate_data, by = c("FID", "IID"), all.x = TRUE)
  
  # Remove rows with missing essential variables
  essential_vars <- c("Outcome", "Add_PRS", "Int_PRS", "Cov_PRS", "Covariate_Pheno")
  complete_cases <- stats::complete.cases(merged_data[, essential_vars])
  merged_data <- merged_data[complete_cases, ]
  
  # Check if there is sufficient data for analysis
  if(nrow(merged_data) < 10) {
    stop("Error: Insufficient complete observations for analysis (n < 10)")
  }
  
  if(verbose) {
    cat(sprintf("Final dataset contains %d observations with %d variables.\n", 
                nrow(merged_data), ncol(merged_data)))
    cat(sprintf("Outcome summary after merging: Mean = %.3f, SD = %.3f\n",
                mean(merged_data$Outcome, na.rm = TRUE),
                stats::sd(merged_data$Outcome, na.rm = TRUE)))
  }
  
  # Identify confounder variables (if any)
  confounder_vars <- colnames(merged_data)[grepl("^Conf_", colnames(merged_data))]
  
  # Build formula string
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS"
  
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }

  model_formula <- stats::as.formula(formula_str)

  if(verbose) {
    cat("Fitting linear regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")
  }

  # Fit linear regression model
  tryCatch({
    model <- stats::lm(model_formula, data = merged_data)
    model_summary <- summary(model)
    
    # Extract model diagnostics
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared
    residual_se <- model_summary$sigma
    f_statistic <- model_summary$fstatistic
    
    # Calculate AIC and BIC
    model_aic <- stats::AIC(model)
    model_bic <- stats::BIC(model)

    return(list(
      model = model,
      model_summary = model_summary,
      merged_data = merged_data,
      add_prs = add_prs,
      int_prs = int_prs,
      cov_prs = cov_prs,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      residual_se = residual_se,
      f_statistic = f_statistic,
      aic = model_aic,
      bic = model_bic,
      scaled = scale_prs,
      temp_files_saved = save_temp_files,
      temp_dir = temp_dir,
      temp_files = temp_files
    ))
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    
    return(list(
      error = e$message,
      merged_data = merged_data,
      add_prs = add_prs,
      int_prs = int_prs,
      cov_prs = cov_prs,
      formula = model_formula,
      temp_dir = temp_dir,
      temp_files = temp_files
    ))
  })
}