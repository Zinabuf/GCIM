#' Perform GWAS for quantitative covariates.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param qp_dis_cov File path for covariate data in the discovery dataset.
#' @param temp_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
q_gwas <- function(plink_path, dis_snp, qp_dis_cov) {
  temp_dir <- tempdir()
system(paste(
  plink_path, "--bfile", dis_snp,
  "--pheno", qp_dis_cov, "--pheno-col-nums", 3,  # Specify that column 3 contains phenotype
  "--glm",
  "--covar", qp_dis_cov, "--covar-col-nums", "4-19",  # Specify covariates
  "--allow-no-sex",
  "--covar-variance-standardize",
  "--out", file.path(temp_dir, "q_gwas"), "&> gwas.log"
))

  covariate_result <- read.table(file.path(temp_dir, "q_gwas.PHENO1.glm.linear"),
                                 header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  covariate_additive <- covariate_result[covariate_result$TEST == "ADD", ]
  write.table(
    covariate_additive[, c("ID", "A1", "BETA")],
    file.path(temp_dir, "covadd.txt"),
    row.names = FALSE, sep = " ", quote = FALSE
  )
}