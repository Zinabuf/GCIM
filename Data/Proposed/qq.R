# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "<plink/path>/plink2"

# For quantitative traits, use corresponding functions
# We conducted a GWEIS of the outcome variable to generate both additive and interaction PRS.
GWEIS <- GWEIS_quantitative(plink_path, "DummyData", "Qphen_disc.txt", "Qexp_dis_cov.txt")
# We conducted a GWAS of the exposure variable to construct an exposure PRS.
GWAS <- GWAS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qcov_disc.txt")

# Extract and compute PRS
# Extract summary statistics, such as extracting the additive and interaction components from the GWEIS for the outcome variable. 
add <- GWEIS[c("ID", "A1", "ADD_BETA")]
int <- GWEIS[c("ID", "A1", "INTERACTION_BETA")]
# Extracting the additive component from the exposure GWAS
add_exp <- GWAS[c("ID", "A1", "BETA")]

add_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
int_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = int)
#Note: If the exposure variable is quantitative, use it as is. If the exposure variable is binary, apply the `PRS_binary` function to generate the object `exp_prs` specified as described below. 
exp_prs <- PRS_quantitative(plink_path, "DummyData", summary_input = add_exp)

# Run GCIM analysis with automatic saving and scaling
# This model specification corresponds to Model 4, as presented in the manuscript.
result1 <- gcim_q0("Qphen_tar.txt", "Qexp_tar_cov.txt", 
                 Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs)
#This model additionally incorporates the main effect of the PRS for the exposure variable. Details of this specification are provided in the model equation section of the Results. 
result2 <- gcim_q1("Qphen_tar.txt", "Qexp_tar_cov.txt", 
                 Add_PRS = add_prs, Int_PRS = int_prs, Cov_PRS = exp_prs)
# Access results
 print(result1$model_summary)
###
print(result2$model_summary)