#' Perform GWEIS for quantitative outcome and quantitative exposure variables.
#'
#' This function performs genome-wide interaction studies (GWEIS) and processes the results
#' to generate files for downstream analysis.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param qp_dis_phen File path for the phenotype data in the discovery dataset.
#' @param qp_dis_cov File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
qqp_gweis <- function(plink_path, dis_snp, qp_dis_phen, qp_dis_cov, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Run PLINK for GWEIS
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", qp_dis_phen,
    "--covar", qp_dis_cov, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "qqp_gweis"), "&> plink.log"
  ))

  # Process GWEIS results
  gweis_result <- read.table(file.path(output_dir, "qqp_gweis.PHENO1.glm.linear"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  # Save additive and interaction results
  additive_outcome <- gweis_result[gweis_result$TEST == "ADD", ]
  write.table(
    additive_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "phenadd_qqp.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_qqp.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )
}
#' Perform GWAS for covariates.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param qp_dis_cov File path for covariate data in the discovery dataset.
#' @param output_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
qp_gwas <- function(plink_path, dis_snp, qp_dis_cov, output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  system(paste(
    plink_path, "--bfile", dis_snp,
    "--pheno", qp_dis_cov, "--glm",
    "--covar-col-nums 4-19 --allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "qp_gwas"), "&> plink.log"
  ))

  covariate_result <- read.table(file.path(output_dir, "qp_gwas.PHENO1.glm.linear"),
                                 header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  covariate_additive <- covariate_result[covariate_result$TEST == "ADD", ]
  write.table(
    covariate_additive[, c("ID", "A1", "BETA")],
    file.path(output_dir, "covadd_qp.txt"),
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
qqp_prs <- function(plink_path, tar_snp, output_dir) {
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  # Compute PRS for additive, interaction, and covariate scores
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_qqp.txt"),
               "--out", file.path(output_dir, "add_qqp")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_qqp.txt"),
               "--out", file.path(output_dir, "int_qqp")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_qp.txt"),
               "--out", file.path(output_dir, "covadd_qp")))

  # Scale PRS values
  prs_add <- scale(read.table(file.path(output_dir, "add_qqp.sscore"))$SCORE1_AVG)
  prs_int <- scale(read.table(file.path(output_dir, "int_qqp.sscore"))$SCORE1_AVG)
  prs_cov <- scale(read.table(file.path(output_dir, "covadd_qp.sscore"))$SCORE1_AVG)

  # Save scaled PRS values
  write.table(prs_add, file.path(output_dir, "add_qqp_scaled.txt"), row.names = FALSE, col.names = FALSE)
  write.table(prs_int, file.path(output_dir, "int_qqp_scaled.txt"), row.names = FALSE, col.names = FALSE)
  write.table(prs_cov, file.path(output_dir, "covadd_qp_scaled.txt"), row.names = FALSE, col.names = FALSE)

  return(list(Additive = prs_add, Interaction = prs_int, Covariate = prs_cov))
}
#' Perform Regression Analysis for GCIM.
#'
#' @param qp_tar_phen File path for the target phenotype data.
#' @param qp_tar_cov File path for the target covariate data.
#' @param prs_add_scaled Scaled additive PRS values.
#' @param prs_int_scaled Scaled interaction PRS values.
#' @param prs_cov_scaled Scaled covariate PRS values.
#' @param confounders Data frame of additional confounders.
#' @return Summary of the regression model.
#' @export
gcim_qqp <- function(qp_tar_phen, qp_tar_cov, prs_add_scaled, prs_int_scaled, prs_cov_scaled, confounders) {
  # Prepare regression data
  regression_data <- data.frame(
    Outcome = read.table(qp_tar_phen)$V3,
    Additive_PRS = prs_add_scaled,
    Interaction_PRS = prs_int_scaled,
    Cov_PRS = prs_cov_scaled,
    Covariate_Pheno = read.table(qp_tar_cov)$V3,
    Confounders = confounders
  )

  # Fit the regression model
  model <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Cov_PRS +
                 Interaction_PRS:Cov_PRS + Confounders,
               family = "binomial", data = regression_data)
  
  # Return model summary
  return(summary(model))
}
