# R/gcim_quantitative.R

#' Perform regression analysis for genetic causality inference model(GCIM) with quantitative outcome
#'
#' This function performs linear regression analysis for GCIM 
#' with quantitative outcomes using without the main effects of Polygenic Risk Scores (PRS) of the exposure.
#' It can read PRS files from temporary directories or use provided PRS data objects.
#'
#' @param qp_tar_phen File path for the target phenotype data (FID, IID, Outcome format)
#' @param qp_tar_cov File path for the target covariate data (FID, IID, Covariate, Confounders format)
#' @param Add_PRS data frame for additive PRS values or NULL to read from temp files
#' @param Int_PRS data frame for interaction PRS values or NULL to read from temp files
#' @param Cov_PRS data frame for covariate PRS values or NULL to read from temp files
#' @param temp_dir Directory path containing PRS files (default: tempdir())
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#' @param scale_prs Logical, whether to scale PRS values (default: TRUE)
#' @param save_temp_files Logical, whether to save temporary PRS files (default: TRUE)
#' @return List containing model summary and diagnostic information
#' @export
#' @examples
#' \dontrun{
#' # Method 1: Using PRS_quantitative outputs directly
#' q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)  # Additive PRS
#' r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS  
#' p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)  # Covariate PRS
#' 
#' result <- gcim_q("Qphe_target.txt", "Qexp_target.txt", 
#'                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
#' 
#' }
gcim_q0 <- function(qp_tar_phen, qp_tar_cov, Add_PRS = NULL, Int_PRS = NULL, Cov_PRS = NULL,
                   temp_dir = tempdir(), verbose = TRUE, scale_prs = TRUE, save_temp_files = TRUE) {
  
   cat("Loading data files for quantitative outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- utils::read.table(qp_tar_phen, header = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is numeric
  outcome_data$Outcome <- as.numeric(outcome_data$Outcome)
  
  # Check for valid quantitative outcome
  if(all(is.na(outcome_data$Outcome))) {
    stop("Error: All outcome values are missing (NA)")
  }
  
    cat(sprintf("Outcome summary: Mean = %.3f, SD = %.3f, Range = [%.3f, %.3f]\n",
                mean(outcome_data$Outcome, na.rm = TRUE),
                stats::sd(outcome_data$Outcome, na.rm = TRUE),
                min(outcome_data$Outcome, na.rm = TRUE),
                max(outcome_data$Outcome, na.rm = TRUE)))
  
  # Function to read PRS file with priority for provided data objects
  read_prs_file <- function(prs_type, provided_data = NULL) {
  cat("Loading", prs_type, "PRS data...\n")
    
    # First, try to use provided data object if available
    if(!is.null(provided_data)) {
      if(is.data.frame(provided_data)) {
        # Handle data frame input
        if(ncol(provided_data) >= 3) {
          prs_data <- provided_data[, c(1, 2, 3)]
          colnames(prs_data) <- c("FID", "IID", "PRS")
          prs_data$PRS <- as.numeric(prs_data$PRS)
          prs_data <- prs_data[!is.na(prs_data$PRS), ]
          
          if(nrow(prs_data) > 0) {
             cat("Using provided data object with", nrow(prs_data), "samples.\n")
            return(prs_data)
          }
        }
      } else if(is.character(provided_data) && length(provided_data) == 1 && file.exists(provided_data)) {
        # Handle file path input
        tryCatch({
          prs_data <- utils::read.table(provided_data, header = TRUE, stringsAsFactors = FALSE)
          if(ncol(prs_data) >= 3) {
            prs_data <- prs_data[, c(1, 2, 3)]
            colnames(prs_data) <- c("FID", "IID", "PRS")
            prs_data$PRS <- as.numeric(prs_data$PRS)
            prs_data <- prs_data[!is.na(prs_data$PRS), ]
            
            if(nrow(prs_data) > 0) {
              if(verbose) cat("Using provided file path with", nrow(prs_data), "samples.\n")
              return(prs_data)
            }
          }
        }, error = function(e) {
          cat("Error reading provided file:", e$message, "\n")
        })
      }
    }
        
    stop("No valid ", prs_type, " PRS data found. Please provide valid data object or ensure files exist in ", temp_dir)
  }
  
  # Read PRS data with priority for provided data objects
  add_prs <- read_prs_file("add", Add_PRS)
  int_prs <- read_prs_file("int", Int_PRS) 
  cov_prs <- read_prs_file("cov", Cov_PRS)
  
  # Rename columns for clarity
  colnames(add_prs) <- c("FID", "IID", "Add_PRS")
  colnames(int_prs) <- c("FID", "IID", "Int_PRS")
  colnames(cov_prs) <- c("FID", "IID", "Cov_PRS")
  
  # Scale the PRS scores if requested
  if(scale_prs) {
     cat("Scaling PRS values (mean=0, sd=1)...\n")
    add_prs$Add_PRS <- as.numeric(scale(add_prs$Add_PRS)[,1])
    int_prs$Int_PRS <- as.numeric(scale(int_prs$Int_PRS)[,1])
    cov_prs$Cov_PRS <- as.numeric(scale(cov_prs$Cov_PRS)[,1])
  }
  
    cat("PRS data summary:\n")
    cat("- Additive PRS: ", nrow(add_prs), " individuals\n")
    cat("- Interaction PRS: ", nrow(int_prs), " individuals\n") 
    cat("- Covariate PRS: ", nrow(cov_prs), " individuals\n")

  
  # Save copies to standardized temporary files if requested
  temp_files <- NULL
  if(save_temp_files) {
    cat("Saving PRS copies to standardized temporary files...\n")
    
    # Define standardized temporary file paths  
    add_prs_file <- file.path(temp_dir, "add_prs_standardized.txt")
    int_prs_file <- file.path(temp_dir, "int_prs_standardized.txt")
    cov_prs_file <- file.path(temp_dir, "cov_prs_standardized.txt")
    
    # Save copies with headers for easier reading
    utils::write.table(add_prs, add_prs_file, quote = FALSE, row.names = FALSE, sep = "\t")
    utils::write.table(int_prs, int_prs_file, quote = FALSE, row.names = FALSE, sep = "\t")
    utils::write.table(cov_prs, cov_prs_file, quote = FALSE, row.names = FALSE, sep = "\t")
    
    # Store file paths for return object
    temp_files <- c(add_prs_file, int_prs_file, cov_prs_file)
    names(temp_files) <- c("add_prs", "int_prs", "cov_prs")
  }
  
  # Load covariate data
  covariate_data <- utils::read.table(qp_tar_cov, header = FALSE)
  
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
   cat("Merging data files by FID and IID...\n")
  
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
  if(nrow(merged_data) < 1) {
    stop("Error: Insufficient complete observations for analysis (n < 1)")
  }
  
  
    cat(sprintf("Final dataset contains %d observations with %d variables.\n", 
                nrow(merged_data), ncol(merged_data)))
    cat(sprintf("Outcome summary after merging: Mean = %.3f, SD = %.3f\n",
                mean(merged_data$Outcome, na.rm = TRUE),
                stats::sd(merged_data$Outcome, na.rm = TRUE)))
  
  
  # Identify confounder variables (if any)
  confounder_vars <- colnames(merged_data)[grepl("^Conf_", colnames(merged_data))]
  
  # Build formula string
  formula_str <- "Outcome ~ Add_PRS + Int_PRS + Covariate_Pheno + Int_PRS:Cov_PRS"
  
  if (length(confounder_vars) > 0) {
    formula_str <- paste(formula_str, "+", paste(confounder_vars, collapse = " + "))
  }

  model_formula <- stats::as.formula(formula_str)


    cat("Fitting linear regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")

  # Fit linear regression model
  tryCatch({
    model <- stats::lm(model_formula, data = merged_data)
    model_summary <- summary(model)
    
    # Extract basic model information
    r_squared <- model_summary$r.squared
    adj_r_squared <- model_summary$adj.r.squared
    residual_se <- model_summary$sigma
    f_statistic <- model_summary$fstatistic
    
    return(list(
      model = model,
      model_summary = model_summary,
      merged_data = merged_data,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      r_squared = r_squared,
      adj_r_squared = adj_r_squared,
      residual_se = residual_se,
      f_statistic = f_statistic,
      temp_files_saved = save_temp_files,
      temp_dir = temp_dir,
      temp_files = temp_files
    ))
  }, error = function(e) {
    warning("Model fitting failed: ", e$message)
    
    return(list(
      error = e$message,
      merged_data = merged_data,
      formula = model_formula,
      temp_dir = temp_dir,
      temp_files = temp_files
    ))
  })
}