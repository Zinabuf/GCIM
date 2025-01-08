#' Perform GWEIS, GWAS, PRS computation, and regression analysis for the reverse direction of causation in binary outcome and binary exposure variable.
#'
#' This function performs genome-wide interaction studies (GWEIS), genome-wide association studies (GWAS),
#' polygenic risk score (PRS) computation, and regression analysis to determine causal directions
#' for binary outcomes and binary exposure variables in the reverse direction.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset (common prefix for .fam, .bim, .bed files).
#' @param tar_snp Prefix for binary files for the target dataset (common prefix for .fam, .bim, .bed files).
#' @param br_dis_phen File path for the phenotype data of the binary exposure in the discovery dataset for the reverse direction.
#' @param br_tar_phen File path for the phenotype data of the binary exposure in the target dataset for the reverse direction.
#' @param br_dis_cov File path for phenotype data in the discovery dataset for binary exposure as an outcome in the reverse direction.
#' @param br_tar_cov File path for phenotype data in the target dataset for binary exposure as an outcome in the reverse direction.
#' @param output_dir Directory for saving output files.
#' @param confounders Covariates to include as confounders in the regression model.
#'
#' @return Summary of the regression model.
#' @export
bbr_gweis <- function(plink_path, dis_snp, tar_snp, br_dis_phen, br_tar_phen, br_dis_cov, br_tar_cov, output_dir, confounders) {
  # Ensure output directory exists
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  }
  # Step 1: Perform GWEIS
  perform_gweis <- function() {
    system(paste(
      plink_path, "--bfile", dis_snp,
      "--glm interaction --pheno", br_dis_cov,
      "--covar", br_dis_phen, "--parameters 1,2,3,4-19",
      "--allow-no-sex --covar-variance-standardize",
      "--out", file.path(output_dir, "bbr_gweis"), "&> plink.log"
    ))
    
    gweis_result <- read.table(file.path(output_dir, "bbr_gweis.PHENO1.glm.logistic.hybrid"),
                               header = TRUE, stringsAsFactors = FALSE, comment.char = "")
    additive_cov <- gweis_result[gweis_result$TEST == "ADD", ]
    additive_cov$BETA <- log(as.numeric(as.character(additive_cov$OR)))
    write.table(
      additive_cov[, c("ID", "A1", "BETA")],
      file.path(output_dir, "covadd_bbr.txt"),
      row.names = FALSE, sep = " ", quote = FALSE
    )
    
    interaction_cov <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
    interaction_cov$BETA <- log(as.numeric(as.character(interaction_cov$OR)))
    write.table(
      interaction_cov[, c("ID", "A1", "BETA")],
      file.path(output_dir, "int_bbr.txt"),
      row.names = FALSE, sep = " ", quote = FALSE
    )
  }
  
  # Step 2: Perform GWAS
  perform_gwas <- function() {
    system(paste(
      plink_path, "--bfile", dis_snp,
      "--pheno", br_dis_phen, "--glm",
      "--covar-col-nums 4-19 --allow-no-sex --covar-variance-standardize",
      "--out", file.path(output_dir, "br_gwas"), "&> plink.log"
    ))
    
    phen_gwas <- read.table(file.path(output_dir, "br_gwas.PHENO1.glm.logistic.hybrid"),
                            header = TRUE, stringsAsFactors = FALSE, comment.char = "")
    phen_additive <- phen_gwas[phen_gwas$TEST == "ADD", ]
    phen_additive$BETA <- log(as.numeric(as.character(phen_additive$OR)))
    write.table(
      phen_additive[, c("ID", "A1", "BETA")],
      file.path(output_dir, "phenadd_br.txt"),
      row.names = FALSE, sep = " ", quote = FALSE
    )
  }
  
  # Step 3: Compute PRS
  compute_prs <- function() {
    system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "covadd_bbr.txt"), "--out", file.path(output_dir, "covadd_bbr")))
    system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "int_bbr.txt"), "--out", file.path(output_dir, "int_bbr")))
    system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(output_dir, "phenadd_br.txt"), "--out", file.path(output_dir, "phenadd_qr")))
  #scale PRS values
    prs_cov_scaled <- scale(read.table(file.path(output_dir, "covadd_bbr.sscore"))$SCORE1_AVG)
    prs_int_scaled <- scale(read.table(file.path(output_dir, "int_bbr.sscore"))$SCORE1_AVG)
    prs_phen_scaled <- scale(read.table(file.path(output_dir, "phenadd_qr.sscore"))$SCORE1_AVG)
   } 
# Step 4: Perform regression analysis
   perform_regression <- function(br_tar_cov, prs_cov_scaled, prs_int_scaled, prs_phen_scaled, br_tar_phen, confounders) {
  # Read the file and extract the relevant column
  outcome_data <- read.table(br_tar_cov, header = TRUE)
  covariate_data <- read.table(br_tar_phen, header = TRUE)

  # Create the regression data frame
  regression_data <- data.frame(
    Outcome = outcome_data$V3,  # Assuming V3 is the correct column
    Additive_PRS = prs_cov_scaled,
    Interaction_PRS = prs_int_scaled,
    Cov_PRS = prs_phen_scaled,
    Covariate_Pheno = covariate_data$V3  # Assuming V3 is the correct column
  )
  
  # Add confounders if applicable
  if (!is.null(confounders)) {
    regression_data <- cbind(regression_data, confounders)
  }
  
  # Perform logistic regression
  gcim_bbr <- glm(
    Outcome ~ Additive_PRS + Interaction_PRS + Covariate_Pheno +
      Interaction_PRS:Cov_PRS + confounders,
    family = "binomial", data = regression_data
  )
  
  return(summary(gcim_bbr))
}
