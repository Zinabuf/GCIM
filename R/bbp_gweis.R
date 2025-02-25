#' Perform GWEIS for binary outcome and exposure variables.
#'
#' This function performs genome-wide interaction studies (GWEIS) and processes the results
#' to generate files for downstream analysis.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param bp_dis_phen File path for the phenotype data in the discovery dataset.
#' @param bp_dis_cov File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
bbp_gweis <- function(plink_path, dis_snp, bp_dis_phen, bp_dis_cov, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Run PLINK for GWEIS
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", bp_dis_phen,
    "--covar", bp_dis_cov, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "bbp_gweis"), "&> gweis.log"
  ))

  # Process GWEIS results
  gweis_result <- read.table(file.path(output_dir, "bbp_gweis.PHENO1.glm.logistic.hybrid"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  # Save additive and interaction results
  additive_outcome <- gweis_result[gweis_result$TEST == "ADD", ]
  additive_outcome$BETA <- log(as.numeric(as.character(additive_outcome$OR)))
  write.table(
    additive_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "phenadd_bbp.txt"),
    row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  interaction_outcome$BETA <- log(as.numeric(as.character(interaction_outcome$OR)))
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_bbp.txt"),
    row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE
  )
}
#' Perform GWAS for covariates.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param bp_dis_cov File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
bp_gwas <- function(plink_path, dis_snp, bp_dis_cov, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  system(paste(
    plink_path, "--bfile", dis_snp,
    "--pheno", bp_dis_cov, "--glm",
    "--covar-col-nums 4-19 --allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "bp_gwas"), "&> gwas.log"
  ))

  covariate_result <- read.table(file.path(output_dir, "bp_gwas.PHENO1.glm.logistic.hybrid"),
                                 header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  covariate_additive <- covariate_result[covariate_result$TEST == "ADD", ]
  covariate_additive$BETA <- log(as.numeric(as.character(covariate_additive$OR)))
  write.table(
    covariate_additive[, c("ID", "A1", "BETA")],
    file.path(output_dir, "covadd_bp.txt"),
    row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE
  )
}
#' Compute Polygenic Risk Scores (PRS).
#'
#' @param plink_path Path to the PLINK executable application.
#' @param tar_snp Prefix for binary files for the target dataset.
#' @param output_dir Directory to save output files.
#' @return A list containing scaled PRS values for additive, interaction, and covariate scores.
#' @export
bbp_prs <- function(plink_path, tar_snp, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Compute PRS for additive, interaction, and covariate scores
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_bbp.txt"),
               "--out", file.path(output_dir, "add_bbp"), "&> add_bbp.log"))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_bbp.txt"),
               "--out", file.path(output_dir, "int_bbp"), "&> int_bbp.log"))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_bp.txt"),
               "--out", file.path(output_dir, "covadd_bp"), "&> covadd_bbp.log"))

  # Check if PRS files exist before reading
  prs_files <- list(
    prs_add = file.path(output_dir, "add_bbp.sscore"),
    prs_int = file.path(output_dir, "int_bbp.sscore"),
    prs_cov = file.path(output_dir, "covadd_bp.sscore")
  )
  prs_values <- list()
  for (name in names(prs_files)) {
    if (file.exists(prs_files[[name]])) {
      prs_data <- read.table(prs_files[[name]], header = FALSE)
      prs_data[, 5] <- scale(prs_data[, 5])  # Scale only column 5
      prs_extracted <- prs_data[, c(1, 2, 5)]  # Extract columns 1, 2, and 5
      write.table(prs_extracted, file.path(output_dir, paste0(name, "_scaled.txt")), 
                  row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
       colnames(prs_extracted) = c("FID", "IID", "PRS")
      prs_values[[name]] <- prs_extracted
    } else {
      warning(paste("Warning: PRS file missing -", prs_files[[name]]))
      prs_values[[name]] <- NULL
    }
  }
  return(prs_values)
}
#' Perform Regression Analysis for GCIM.
#'
#' @param bp_tar_phen File path for the target phenotype data.
#' @param bp_tar_cov File path for the target covariate data.
#' @param Additive Scaled additive PRS values.
#' @param Interaction Scaled interaction PRS values.
#' @param Covariate Scaled covariate PRS values.
#' @param Confounders Data frame of additional confounders.
#' @return Summary of the regression model.
#' @export
gcim_bbp <- function(bp_tar_phen, bp_tar_cov, Add_PRS, Int_PRS, Cov_PRS, confounders ) {
  # Load phenotype and covariate data
  outcome_bp_data <- read.table("bp_tar_phen", header = TRUE, stringsAsFactors = FALSE)
  prs_add <- read.table("prs_add_scaled.txt", header = TRUE, stringsAsFactors = FALSE)
  prs_int <- read.table("prs_int_scaled.txt", header = TRUE, stringsAsFactors = FALSE)
  prs_cov <- read.table("prs_cov_scaled.txt", header = TRUE, stringsAsFactors = FALSE)
  covariate_bp_data <- read.table("bp_tar_cov", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

  # Prepare regression data
  regression_data <- data.frame(
    Outcome = outcome_bp_data[, 3],  # Assuming 3rd column is the outcome
    Add_PRS = prs_add[, 3],   # Adjust column index if necessary
    Int_PRS = prs_int[, 3],
    Cov_PRS = prs_cov[, 3],
    Covariate_Pheno = covariate_bp_data[, 3] # Adjust column index as needed
  )

  # Add confounders dynamically based on the input `conf_`
  confounders <- covariate_bp_data[, 4:19]
  colnames(confounders) <- paste0("Conf_", seq_along(confounders))  # Rename confounders

  # Combine confounders with regression data
  regression_data <- cbind(regression_data, confounders)

  # Fit the regression model using all variables
  model_formula <- as.formula(paste("Outcome ~ Add_PRS + Int_PRS + Cov_PRS + Int_PRS:Cov_PRS +", paste(names(confounders), collapse = " + ")))
  model <- glm(model_formula, family = binomial(), data = regression_data)

  # Return model summary
  return(summary(model))
}