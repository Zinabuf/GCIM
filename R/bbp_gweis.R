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
    "--out", file.path(output_dir, "bbp_gweis"), "&> plink.log"
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
    row.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  interaction_outcome$BETA <- log(as.numeric(as.character(interaction_outcome$OR)))
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_bbp.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
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
    "--out", file.path(output_dir, "bp_gwas"), "&> plink.log"
  ))

  covariate_result <- read.table(file.path(output_dir, "bp_gwas.PHENO1.glm.logistic.hybrid"),
                                 header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  covariate_additive <- covariate_result[covariate_result$TEST == "ADD", ]
  covariate_additive$BETA <- log(as.numeric(as.character(covariate_additive$OR)))
  write.table(
    covariate_additive[, c("ID", "A1", "BETA")],
    file.path(output_dir, "covadd_bp.txt"),
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
bbp_prs <- function(plink_path, tar_snp, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Compute PRS for additive, interaction, and covariate scores
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_bbp.txt"),
               "--out", file.path(output_dir, "add_bbp")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_bbp.txt"),
               "--out", file.path(output_dir, "int_bbp")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_bp.txt"),
               "--out", file.path(output_dir, "covadd_bp")))

  # Scale PRS values
  prs_add <- scale(read.table(file.path(output_dir, "add_bbp.sscore"))$SCORE1_AVG)
  prs_int <- scale(read.table(file.path(output_dir, "int_bbp.sscore"))$SCORE1_AVG)
  prs_cov <- scale(read.table(file.path(output_dir, "covadd_bp.sscore"))$SCORE1_AVG)

  return(list(Additive = prs_add, Interaction = prs_int, Covariate = prs_cov))
}
#' Perform Regression Analysis for GCIM.
#'
#' @param bp_tar_phen File path for the target phenotype data.
#' @param bp_tar_cov File path for the target covariate data.
#' @param prs_add Scaled additive PRS values.
#' @param prs_int Scaled interaction PRS values.
#' @param prs_cov Scaled covariate PRS values.
#' @param confounders Data frame of additional confounders.
#' @return Summary of the regression model.
#' @export
gcim_bbp <- function(bp_tar_phen, bp_tar_cov, Additive, Interaction, Covariate, confounders) {
  # Prepare regression data
  regression_data <- data.frame(
    Outcome = read.table(bp_tar_phen)$V3,
    Additive_PRS = prs_add,
    Interaction_PRS = prs_int,
    Cov_PRS = prs_cov,
    Covariate_Pheno = read.table(bp_tar_cov)$V3,
    Confounders = confounders
  )

  # Fit the regression model
  model <- glm(Outcome ~ Additive + Interaction + Covariate_Pheno +
                 Interaction_PRS:Covariate,
               family = "binomial", data = regression_data)
  
  # Return model summary
  return(summary(model))
}