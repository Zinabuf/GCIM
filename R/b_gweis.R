#' Perform GWEIS for binary outcome variables.
#'
#' This function performs genome-wide-by-environemnt interaction studies (GWEIS) and processes the results
#' to generate files for downstream analysis.
#'
#' @param plink_path Path to the PLINK executable application.
#' @param dis_snp Prefix for binary files for the discovery dataset.
#' @param bp_dis_phen File path for the phenotype data in the discovery dataset.
#' @param bp_dis_cov File path for covariate data in the discovery dataset.
#' @param temp_dir Directory to save output files.
#' @return None. Results are saved to files.
#' @export
b_gweis <- function(plink_path, dis_snp, bp_dis_phen, bp_dis_cov) {
  temp_dir <- tempdir()
  # Run PLINK for GWEIS
  system(paste(
    plink_path, "--bfile", dis_snp,
    "--glm interaction --pheno", bp_dis_phen,
    "--covar", bp_dis_cov, "--parameters 1,2,3,4-19",
    "--allow-no-sex --covar-variance-standardize",
    "--out", file.path(temp_dir, "b_gweis"), "&> gweis.log"
  ))

  # Process GWEIS results
  gweis_result <- read.table(file.path(temp_dir, "b_gweis.PHENO1.glm.logistic.hybrid"),
                             header = TRUE, stringsAsFactors = FALSE, comment.char = "")

  # Save additive and interaction results
  additive_outcome <- gweis_result[gweis_result$TEST == "ADD", ]
  additive_outcome$BETA <- log(as.numeric(as.character(additive_outcome$OR)))
  write.table(
    additive_outcome[, c("ID", "A1", "BETA")],
    file.path(temp_dir, "phenadd.txt"),
    row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE
  )

  interaction_outcome <- gweis_result[gweis_result$TEST == "ADDxCOVAR1", ]
  interaction_outcome$BETA <- log(as.numeric(as.character(interaction_outcome$OR)))
  write.table(
    interaction_outcome[, c("ID", "A1", "BETA")],
    file.path(temp_dir, "int.txt"),
    row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE
  )
}
