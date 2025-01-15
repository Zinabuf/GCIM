#' Perform GWEIS for quantitative outcome and binary exposure variables for the reverse directions.
#'
#' This function performs genome-wide interaction studies (GWEIS) and processes the results
#' to generate files for downstream analysis.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param qr_dis_phen File path for the phenotype data in the discovery dataset.
#' @param br_dis_cov File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
qbr_gweis <- function(plink_path, dis_snp, br_dis_cov, qr_dis_phen, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Run PLINK for GWEIS
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", br_dis_cov,
    "--covar", qr_dis_phen, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "qbr_gweis"), "&> plink.log"
  ))

  # Process GWEIS results
  gweis_result <- read.table(file.path(output_dir, "qbr_gweis.PHENO1.glm.logistic.hybrid"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  # Save additive and interaction results
  additive_outcome <- gweis_result[gweis_result$TEST == "ADD", ]
  additive_outcome$BETA <- log(as.numeric(as.character(additive_outcome$OR)))
  write.table(
    additive_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "phenadd_qbr.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  interaction_outcome$BETA <- log(as.numeric(as.character(interaction_outcome$OR)))
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_qbr.txt"),
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
qbr_prs <- function(plink_path, tar_snp, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Compute PRS for additive, interaction, and covariate scores
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_qbr.txt"),
               "--out", file.path(output_dir, "add_qbr")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_qbr.txt"),
               "--out", file.path(output_dir, "int_qbr")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_qr.txt"),
               "--out", file.path(output_dir, "covadd_qr")))

  # Scale PRS values
  prs_add <- scale(read.table(file.path(output_dir, "add_qbr.sscore"))$SCORE1_AVG)
  prs_int <- scale(read.table(file.path(output_dir, "int_qbr.sscore"))$SCORE1_AVG)
  prs_cov <- scale(read.table(file.path(output_dir, "covadd_qr.sscore"))$SCORE1_AVG)

  # Save scaled PRS values
  write.table(prs_add, file.path(output_dir, "add_qbr_scaled.txt"), row.names = FALSE, col.names = FALSE)
  write.table(prs_int, file.path(output_dir, "int_qbr_scaled.txt"), row.names = FALSE, col.names = FALSE)
  write.table(prs_cov, file.path(output_dir, "covadd_qr_scaled.txt"), row.names = FALSE, col.names = FALSE)

  return(list(Additive = prs_add, Interaction = prs_int, Covariate = prs_cov))
}
#' Perform Regression Analysis for GCIM.
#'
#' @param qr_tar_phen File path for the target phenotype data.
#' @param br_tar_cov File path for the target covariate data.
#' @param prs_add_scaled Scaled additive PRS values.
#' @param prs_int_scaled Scaled interaction PRS values.
#' @param prs_cov_scaled Scaled covariate PRS values.
#' @param confounders Data frame of additional confounders.
#' @return Summary of the regression model.
#' @export
gcim_qbr <- function(br_tar_cov, qr_tar_phen, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders) {
  # Prepare regression data
  regression_data <- data.frame(
    Outcome = read.table(br_tar_cov)$V3,
    Additive_PRS = prs_add_scaled,
    Interaction_PRS = prs_int_scaled,
    Cov_PRS = prs_cov_scaled,
    Covariate_Pheno = read.table(qr_tar_phen)$V3,
    Confounders = confounders
  )

  # Fit the regression model
  model <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                 Interaction_PRS:Cov_PRS + Confounders,
               family = "binomial", data = regression_data)
  
  # Return model summary
  return(summary(model))
}
