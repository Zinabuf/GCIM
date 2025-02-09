#' Perform GWEIS for quantitative outcome and quantitative exposure variables for the reverse direction.
#'
#' This function performs genome-wide interaction studies (GWEIS) and processes the results
#' to generate files for downstream analysis.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param qr_dis_phen File path for the phenotype data in the discovery dataset.
#' @param qr_dis_cov File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
qqr_gweis <- function(plink_path, dis_snp, qr_dis_cov, qr_dis_phen, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Run PLINK for GWEIS
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", qr_dis_cov,
    "--covar", qr_dis_phen, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "qqr_gweis"), "&> plink.log"
  ))

  # Process GWEIS results
  gweis_result <- read.table(file.path(output_dir, "qqr_gweis.PHENO1.glm.linear"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  # Save additive and interaction results
  additive_outcome <- gweis_result[gweis_result$TEST == "ADD", ]
  write.table(
    additive_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "phenadd_qqr.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_qqr.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )
}
#' Perform GWAS for covariates.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param qr_dis_phen File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
qr_gwas <- function(plink_path, dis_snp, qr_dis_phen, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  system(paste(
    plink_path, "--bfile", dis_snp,
    "--pheno", qr_dis_phen, "--glm",
    "--covar-col-nums 4-19 --allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "qr_gwas"), "&> plink.log"
  ))

  covariate_result <- read.table(file.path(output_dir, "qr_gwas.PHENO1.glm.linear"),
                                 header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  covariate_additive <- covariate_result[covariate_result$TEST == "ADD", ]
  write.table(
    covariate_additive[, c("ID", "A1", "BETA")],
    file.path(output_dir, "covadd_qr.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )
}
#' Compute Polygenic Risk Scores (PRS).
#'
#' @param plink_path Path to the PLINK executable application.
#' @param tar_snp Prefix for binary files for the target dataset.
#' @param output_dir Directory to save output files.
#' @return A list containing scaled PRS values for additive, interaction, and covariate scores.
#' @export
qqr_prs <- function(plink_path, tar_snp, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Compute PRS for additive, interaction, and covariate scores
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_qqr.txt"),
               "--out", file.path(output_dir, "add_qqr"), "&> add_qqr.log"))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_qqr.txt"),
               "--out", file.path(output_dir, "int_qqr"), "&> int_qqr.log"))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_qr.txt"),
               "--out", file.path(output_dir, "covadd_qr"), "&> covadd_qr.log"))

  # Check if PRS files exist before reading
  prs_files <- list(
    prs_add = file.path(output_dir, "add_qqr.sscore"),
    prs_int = file.path(output_dir, "int_qqr.sscore"),
    prs_cov = file.path(output_dir, "covadd_qr.sscore")
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
#' @param qr_tar_phen File path for the target phenotype data.
#' @param qr_tar_cov File path for the target covariate data.
#' @param prs_add_scaled Scaled additive PRS values.
#' @param prs_int_scaled Scaled interaction PRS values.
#' @param prs_cov_scaled Scaled covariate PRS values.
#' @param confounders Data frame of additional confounders.
#' @return Summary of the regression model.
#' @export
gcim_qqr <- function(qr_tar_cov, qr_tar_phen, Add_PRS, Int_PRS, Cov_PRS, confounders ) {
  # Load phenotype and covariate data
  covariate_qr_data <- read.table("qr_tar_cov", header = TRUE, stringsAsFactors = FALSE)
  prs_add <- read.table("prs_add_scaled.txt", header = TRUE, stringsAsFactors = FALSE)
  prs_int <- read.table("prs_int_scaled.txt", header = TRUE, stringsAsFactors = FALSE)
  prs_cov <- read.table("prs_cov_scaled.txt", header = TRUE, stringsAsFactors = FALSE)
  outcome_qr_data <- read.table("qr_tar_phen", header = TRUE, stringsAsFactors = FALSE, fill = TRUE)

  # Prepare regression data
  regression_data <- data.frame(
    Outcome = covariate_qr_data[, 3],  # Assuming 3rd column is the outcome
    Add_PRS = prs_add[, 3],   # Adjust column index if necessary
    Int_PRS = prs_int[, 3],
    Cov_PRS = prs_cov[, 3],
    Covariate_Pheno = outcome_qr_data[, 3] # Adjust column index as needed
  )

  # Add confounders dynamically based on the input `conf_`
  confounders <- outcome_qr_data[, 4:19]
  colnames(confounders) <- paste0("Conf_", seq_along(confounders))  # Rename confounders

  # Combine confounders with regression data
  regression_data <- cbind(regression_data, confounders)

  # Fit the regression model using all variables
  model_formula <- as.formula(paste("Outcome ~ Add_PRS + Int_PRS + Cov_PRS + Int_PRS:Cov_PRS +", paste(names(confounders), collapse = " + ")))
  model <- lm(model_formula, data = regression_data)

  # Return model summary
  return(summary(model))
}