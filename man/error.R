bbp_prs <- function(plink_path, tar_snp, output_dir) {
  # Run PLINK commands (assuming these are correct)
  additive_file <- file.path(output_dir, "add_bbp.sscore")
  interaction_file <- file.path(output_dir, "int_bbp.sscore")
  covariate_file <- file.path(output_dir, "covadd_bp.sscore")
  
  # Check if output files exist
  if (!file.exists(additive_file) || !file.exists(interaction_file) || !file.exists(covariate_file)) {
    stop("Error: One or more score files missing. Check PLINK execution.")
  }
  
  # Read scores from PLINK output
  Additive <- read.table(additive_file, header = TRUE)$SCORE
  Interaction <- read.table(interaction_file, header = TRUE)$SCORE
  Covariate <- read.table(covariate_file, header = TRUE)$SCORE
  
  # Return structured output
  return(list(Additive = Additive, Interaction = Interaction, Covariate = Covariate))
}
