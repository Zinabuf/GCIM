# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink/path>/plink2"

# For quantitative traits, use corresponding functions
# We conducted a GWEIS of the outcome variable (treated as the exposure in the proposed causal direction) to generate both additive and interaction PRS.
GWEIS <- GWEIS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qphen_disc_cov.txt")
# We conducted a GWAS of the exposure variable (treated as the outcome in the proposed causal directions) to construct an exposure PRS
GWAS <- GWAS_quantitative(plink_path, "DummyData", "Qphen_disc.txt", "Qcov_disc.txt")

# Extract and compute PRS

add_exp <- GWEIS[c("ID", "A1", "ADD_BETA")]
int_exp <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]

add_out <- GWAS[c("ID", "A1", "BETA")]

add_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add_exp)
int_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = int_exp) 
out_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add_out)

# Run GCIM analysis with automatic saving and scaling
 result1 <- gcim_q0("Qexp_tar.txt", "Qphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs)
 result2 <- gcim_q1("Qexp_tar.txt", "Qphen_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = out_prs)
 # Access results
 print(result1$model_summary)
 print(result2$model_summary)