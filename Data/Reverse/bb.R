# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink/path>/plink2"
# Run GxEprs analysis for binary traits

# We conducted a GWEIS of the outcome variable (treated as the exposure in the proposed causal direction) to generate both additive and interaction PRS.
GWEIS <- GWEIS_binary(plink_path, "DummyData", "Bexp_disc.txt", "Bphen_disc_cov.txt")
# We conducted a GWAS of the exposure variable (treated as the outcome in the proposed causal directions) to construct an exposure PRS
GWAS <- GWAS_binary(plink_path, "DummyData", "Bphen_disc.txt", "Bcov_disc.txt")

# Extract summary statistics
add_exp <- GWEIS[c("ID", "A1", "ADD_BETA")]
int_exp <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]
add_out <- GWAS[c("ID", "A1", "BETA")]

# Compute PRS for each component
add_prs <- PRS_binary(plink_path, "DummyData", summary_input = add_exp)  # Additive PRS
int_prs <- PRS_binary(plink_path, "DummyData", summary_input = int_exp)  # Interaction PRS
out_prs <- PRS_binary(plink_path, "DummyData", summary_input = add_prs)  # Covariate PRS

# Run GCIM analysis with automatic saving and scaling
# A similar model specification is applied as described above for the quantitative analysis.
 result1 <- gcim_b0("Bexp_tar.txt", "Bphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs) 
 result2 <- gcim_b1("Bexp_tar.txt", "Bphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs)
# Access results
 print(result1$model_summary)
 print(result2$model_summary)