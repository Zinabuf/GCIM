# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink/path>/plink2"
# Run GxEprs analysis for binary traits
#  We conducted a GWEIS of the outcome variable to generate both additive and interaction PRS.
GWEIS <- GWEIS_binary(plink_path, "DummyData", "Bphen_disc.txt", "Bexp_disc_cov.txt")
# Note: If the exposure variable is binary, use it as is. If the exposure variable is quantitative, apply the `GWAS_quantitative` function to generate the object `GWAS` specified as described below.
# We conducted a GWAS of the exposure variable to construct an exposure PRS.
GWAS <- GWAS_binary(plink_path, "DummyData", "Bexp_disc.txt", "Bcov_disc.txt")

# Extract summary statistics
# Extracting the additive and interaction components from the GWEIS for the outcome variable. 
add <- GWEIS[c("ID", "A1", "ADD_BETA")]
int <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]
# Extracting the additive component from the exposure GWAS
add_exp <- GWAS[c("ID", "A1", "BETA")]

# Compute PRS for each component
add_prs <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
int_prs <- PRS_binary(plink_path, "DummyData", summary_input = int)  # Interaction PRS
#Note: If the exposure variable is binary, use it as is. If the exposure variable is quantitative, apply the `PRS_quantitative` function to generate the object `exp_prs` specified as described below.
exp_prs <- PRS_binary(plink_path, "DummyData", summary_input = add_exp)  # Covariate PRS

# Run GCIM analysis with automatic saving and scaling
# A similar model specification is applied as described above for the quantitative analysis.
 result1 <- gcim_b0("Bphen_tar.txt", "Bexp_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs) 
 result2 <- gcim_b1("Bphen_tar.txt", "Bexp_tar_cov.txt", 
                  Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs)
# Access results
 print(result1$model_summary)
###
print(result2$model_summary)