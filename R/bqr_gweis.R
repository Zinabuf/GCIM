#' Perform GWEIS, GWAS, PRS computation, and regression analysis for the reverse direction of causations in binary outcome and quantitative exposure variable.
#'
#' This function performs genome-wide interaction studies (GWEIS), genome-wide association studies (GWAS),
#' polygenic risk score (PRS) computation, and regression analysis to determine causal directions
#' for binary outcomes and quatitative exposure variables.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset (common prefix for .fam, .bim, .bed files).
#' @param tar_snp Prefix for binary files for the target dataset (common prefix for .fam, .bim, .bed files).
#' @param br_dis_phen File path for the phenotype data of the binary exposure in the discovery dataset for the reverse direction.
#' @param br_tar_phen File path for the phenotype data of the binary exposure in the target dataset for the reverse direction.
#' @param qr_dis_cov File path for phenotype data in the discovery dataset for binary exposure as an outcome in the reverse direction.
#' @param qr_tar_cov File path for phenotype data in the target dataset for binary exposure as an outcome in the reverse direction.
#' @param output_dir Directory for saving output files.
#' @param confounders Covariates to include as confounders in the regression model.
#'
#' @return Summary of the regression model.
#' @export
bqr_gweis <- function(plink_path, dis_snp, qr_dis_cov, br_dis_phen, output_dir, confounders) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
}
br_gwas <- function(plink_path, dis_snp, br_dis_phen, output_dir, confounders) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
}

  # Step 1: GWEIS for the outcome
bqr_gweis <- function() {
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", br_dis_cov,
    "--covar", br_dis_phen, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "bqr_gweis"), "&> plink.log"
  ))

  # Load GWEIS results and process
  gweis_result <- read.table(file.path(output_dir, "bqr_gweis.PHENO1.glm.linear"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  additive_cov <- gweis_result[gweis_result$TEST == "ADD", ]
  additive_cov$BETA <- log(as.numeric(as.character(additive_cov$OR)))
  write.table(
    additive_cov[, c("ID", "A1", "BETA")],
    file.path(output_dir, "covadd_bqr.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_cov <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  interaction_cov$BETA <- log(as.numeric(as.character(interaction_cov$OR)))
  write.table(
    interaction_cov[, c("ID", "A1", "BETA")],
    file.path(output_dir, "int_bqr.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )
}
  # Step 2: GWAS for covariate
br_gwas <- function() {
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--pheno", br_dis_phen, "--glm",
    "--covar-col-nums 4-19 --allow-no-sex --covar-variance-standardize",
    "--out", file.path(output_dir, "br_gwas"), "&> plink.log"
  ))

  phen_gwas <- read.table(file.path(output_dir, "br_gwas.PHENO1.glm.logistic.hybrid"),
                                 header = TRUE, stringsAsFactors = FALSE, comment.char = "")
  phen_additive <- phen_result[phen_result$TEST == "ADD", ]
  phen_additive$BETA <- log(as.numeric(as.character(phen_additive$OR)))
  write.table(
    phen_additive[, c("ID", "A1", "BETA")],
    file.path(output_dir, "phenadd_br.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )

 }
  # Step 3: Compute PRS
compute_prs <- function() {
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_bqr.txt"), "--out", file.path(output_dir, "covadd_bqr")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_bqr.txt"), "--out", file.path(output_dir, "int_bqr")))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_br.txt"), "--out", file.path(output_dir, "phenadd_br")))

  # Scale PRS values
  prs_cov_scaled <- scale(read.table(file.path(output_dir, "covadd_bqr.sscore"))$SCORE1_AVG)
  prs_int_scaled <- scale(read.table(file.path(output_dir, "int_bqr.sscore"))$SCORE1_AVG)
  prs_phen_scaled <- scale(read.table(file.path(output_dir, "phenadd_br.sscore"))$SCORE1_AVG)
}
  # Step 4: Perform regression analysis
perform_regression <- function() {
  regression_data <- data.frame(
    Outcome = read.table(qr_tar_cov)$V3,
    Additive_PRS = prs_cov_scaled,
    Interaction_PRS = prs_int_scaled,
    Cov_PRS = prs_phen_scaled,
    Covariate_Pheno = read.table(br_tar_phen)$V3,
    Confounders = confounders
  )

  gcim_bqr <- lm(Outcome ~ Additive_PRS + Interaction_PRS + Covariate_Pheno +
                   Interaction_PRS:Cov_PRS + Confounders,
                  data = regression_data)

  return(summary(gcim_bqr))
}
