# R/gcim_binary.R

#' Perform regression analysis for genetic causality inference model(GCIM) with binary outcome
#'
#' This function performs logistic regression analysis for GCIM 
#' with binary outcomes using without the main effects of Polygenic Risk Scores (PRS) of the exposure.
#' It can read PRS files from temporary directories or use provided PRS data objects.
#'
#' @param bp_tar_phen File path for the target phenotype data (FID, IID, Outcome format)
#' @param bp_tar_cov File path for the target covariate data (FID, IID, Covariate, Confounders format)
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
#' # Method 1: Using PRS_binary outputs directly
#' q <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
#' r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS  
#' p <- PRS_binary(plink_path, "DummyData", summary_input = trd)  # Covariate PRS
#' 
#' result <- gcim_b0("Bphe_target.txt", "Bexp_target.txt", 
#'                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
#' 
#' }
gcim_b0 <- function(bp_tar_phen, bp_tar_cov, Add_PRS = NULL, Int_PRS = NULL, Cov_PRS = NULL,
                   temp_dir = tempdir(), verbose = TRUE, scale_prs = TRUE, save_temp_files = TRUE) {
  
   cat("Loading data files for binary outcome analysis...\n")
  
  # Load phenotype data
  outcome_data <- utils::read.table(bp_tar_phen, header = FALSE, stringsAsFactors = FALSE)
  colnames(outcome_data) <- c("FID", "IID", "Outcome")
  
  # Ensure outcome is binary (0/1)
  if(!all(outcome_data$Outcome %in% c(0, 1, NA))) {
    cat("Warning: Outcome variable contains non-binary values. Converting to binary.\n")
  }
  
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
              cat("Using provided file path with", nrow(prs_data), "samples.\n")
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
  covariate_data <- utils::read.table(bp_tar_cov, header = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  
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

  model_formula <- stats::as.formula(formula_str)

   
    cat("Fitting logistic regression model...\n")
    cat("Formula:", deparse(model_formula), "\n")

  # Fit logistic regression model
  tryCatch({
    model <- stats::glm(model_formula, data = merged_data, family = stats::binomial(link = "logit"))
    model_summary <- summary(model)
    
    # Calculate pseudo R-squared (McFadden's R-squared)
    null_deviance <- model$null.deviance
    residual_deviance <- model$deviance
    pseudo_r_squared <- 1 - (residual_deviance / null_deviance)

    return(list(
      model = model,
      model_summary = model_summary,
      merged_data = merged_data,
      formula = model_formula,
      sample_size = nrow(merged_data),
      variables = colnames(merged_data),
      pseudo_r_squared = pseudo_r_squared,
      null_deviance = null_deviance,
      residual_deviance = residual_deviance,
      aic = model$aic,
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