#' Perform GWEIS, GWAS, PRS computation, and regression analysis for the proposed direction of causations in binary outcome and quantitative exposure variable.
#'
#' This function performs genome-wide interaction studies (GWEIS), genome-wide association studies (GWAS),
#' polygenic risk score (PRS) computation, and regression analysis to determine causal directions
#' for binary outcomes and quantitative exposure variables for the proposed direction.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset (common prefix for .fam, .bim, .bed files).
#' @param tar_snp Prefix for binary files for the target dataset (common prefix for .fam, .bim, .bed files).
#' @param bp_dis_phen File path for the phenotype data of the binary outcome in the discovery dataset.
#' @param bp_tar_phen File path for the phenotype data of the binary outcome in the target dataset.
#' @param qp_dis_cov File path for covariate data (exposure and confounders) in the discovery dataset for quatitative exposure.
#' @param qp_tar_cov File path for covariate data (exposure and confounders) in the target dataset for quatitative exposure.
#' @param output_dir Directory for saving output files.
#' @param confounders Covariates to include as confounders in the regression model.
#'
#' @return Summary of the regression model.
#' @export
bqp_gweis <- function(plink_path, dis_snp, bp_dis_phen, qp_dis_cov, output_dir, confounders) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
}
qp_gwas <- function(plink_path, dis_snp, qp_dis_cov, output_dir, confounders) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
}

  # Step 1: GWEIS for the outcome
bbp_gweis <- function() {
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", bp_dis_phen,
    "--covar", qp_dis_cov, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "bqp_gweis"), "&> plink.log"
  ))

  # Load GWEIS results and process
  gweis_result <- read.table(file.path(output_dir, "bqp_gweis.PHENO1.glm.logistic.hybrid"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  additive_outcome <- gweis_result[gweis_result$TEST == "ADD", ]
  additive_outcome$BETA <- log(as.numeric(as.character(additive_outcome$OR)))
  write.table(
    additive_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "phenadd_bqp.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  interaction_outcome$BETA <- log(as.numeric(as.character(interaction_outcome$OR)))
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_bqp.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )
}
  # Step 2: GWAS for covariate
qp_gwas <- function() {
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
  # Step 3: Compute PRS
compute_prs <- function() {
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_bqp.txt"), "--out", file.path(output_dir, "add_bqp")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_bqp.txt"), "--out", file.path(output_dir, "int_bqp")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_qp.txt"), "--out", file.path(output_dir, "covadd_qp")))

  # Scale PRS values
  prs_add_scaled <- scale(read.table(file.path(output_dir, "add_bqp.sscore"))$SCORE1_AVG)
  prs_int_scaled <- scale(read.table(file.path(output_dir, "int_bqp.sscore"))$SCORE1_AVG)
  prs_cov_scaled <- scale(read.table(file.path(output_dir, "covadd_qp.sscore"))$SCORE1_AVG)
}
  # Step 4: Perform regression analysis
perform_regression <- function() {
  regression_data <- data.frame(
    Outcome = read.table(bp_tar_phen)$V3,
    Additive_PRS = prs_add_scaled,
    Interaction_PRS = prs_int_scaled,
    Cov_PRS = prs_cov_scaled,
    Covariate_Pheno = read.table(qp_tar_cov)$V3,
    Confounders = confounders
  )

  gcim_bqp <- glm(Outcome ~ Additive_PRS + Interaction_PRS + Covariate_Pheno +
                   Interaction_PRS:Cov_PRS + Confounders,
                 family = "binomial", data = regression_data)

  return(summary(gcim_bqp))
}
