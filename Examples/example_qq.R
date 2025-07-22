## For Continuous Outcomes
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "/data/alh-admzw/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qcov_disc.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphe_discovery.txt", "Qcov_discovery.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# Run GCIM for continuous outcome
result_qnt <- gcim_q("Qphe_target.txt", "Bexp_target.txt",
  Add_PRS,
  Int_PRS,
  Cov_PRS
)
print(result_qnt$model_summary)
cat("R-squared:", round(result_cont$r_squared, 4), "\n")