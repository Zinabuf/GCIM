#' Compute Polygenic Risk Scores (PRS).
#'
#' @param plink_path Path to the PLINK executable application.
#' @param tar_snp Prefix for binary files for the target dataset.
#' @param temp_dir Directory to save output files.
#' @return A list containing scaled PRS values for additive, interaction, and covariate scores.
#' @export
prs_scores <- function(plink_path, tar_snp) {
  # Create temporary directory if it doesn't exist
  temp_dir <- tempdir()
  # Compute PRS for additive, interaction, and covariate scores
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(temp_dir, "phenadd.txt"),
               "--out", file.path(temp_dir, "add"), "&> add_bbp.log"))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(temp_dir, "int.txt"),
               "--out", file.path(temp_dir, "int"), "&> int_bbp.log"))
  system(paste(plink_path, "--bfile", tar_snp, "--score", file.path(temp_dir, "covadd.txt"),
               "--out", file.path(temp_dir, "covadd"), "&> covadd_bbp.log"))

  # Check if PRS files exist before reading
  prs_files <- list(
    prs_add = file.path(temp_dir, "add.sscore"),
    prs_int = file.path(temp_dir, "int.sscore"),
    prs_cov = file.path(temp_dir, "covadd.sscore")
  )
  prs_values <- list()
  for (name in names(prs_files)) {
    if (file.exists(prs_files[[name]])) {
      prs_data <- read.table(prs_files[[name]], header = FALSE)
      prs_data[, 5] <- scale(prs_data[, 5])  # Scale only column 5
      prs_extracted <- prs_data[, c(1, 2, 5)]  # Extract columns 1, 2, and 5
      write.table(prs_extracted, file.path(temp_dir, paste0(name, "_scaled.txt")), 
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