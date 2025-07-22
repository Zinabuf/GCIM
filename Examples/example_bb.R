## Complete Workflow Example
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "/data/alh-admzw/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bexp_disc.txt", "Bcov_discovery.txt")
b <- GWEIS_binary(plink_path, "DummyData", "Bphe_discovery.txt", "Bcov_discovery.txt")

# Step 2: Extract summary statistics
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

# Step 3: Compute PRS for each component
q <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS
p <- PRS_binary(plink_path, "DummyData", summary_input = trd)  # Covariate PRS

# Step 4: Run GCIM analysis using the convenience function
result <- gcim_b("Bphe_target.txt", "Bexp_target.txt", 
  Add_PRS,
  Int_PRS,
  Cov_PRS,
  verbose = TRUE
)

# Step 5: Examine results
print(result$model_summary)
cat("Sample size:", result$sample_size, "\n")
cat("Cases:", result$cases, "\n") 
cat("Controls:", result$controls, "\n")
cat("Pseudo R-squared:", round(result$pseudo_r2, 4), "\n")

# View odds ratios for key effects
print("Odds Ratios:")
print(round(result$odds_ratios, 3))