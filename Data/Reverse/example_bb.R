## Complete Workflow Example
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "/data/alh-admzw/plink2"
# Step 1: Run GxEprs analysis for binary traits
a <- GWAS_binary(plink_path, "DummyData", "Bphen_disc.txt", "Bcov_disc.txt")
b <- GWEIS_binary(plink_path, "DummyData", "Bexp_disc.txt", "Bphen_disc_cov.txt")

# Step 2: Extract summary statistics
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

# Step 3: Compute PRS for each component
q <- PRS_binary(plink_path, "DummyData", summary_input = add)  # Additive PRS
r <- PRS_binary(plink_path, "DummyData", summary_input = gxe)  # Interaction PRS
p <- PRS_binary(plink_path, "DummyData", summary_input = trd)  # Covariate PRS

# Step 4: Run GCIM analysis with automatic saving and scaling
 result0 <- gcim_b0("Bexp_tar.txt", "Bphen_tar_cov.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
 result1 <- gcim_b1("Bexp_tar.txt", "Bphen_tar_cov.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p) 
# # Step 5: Access results and processed PRS objects
 print(result0$model_summary)
 print(result1$model_summary)
 cat("Sample size:", result0$sample_size, "\n")
 cat("Pseudo R-squared:", result0$pseudo_r_squared, "\n")
 cat("AIC:", result0$aic, "\n")
# 
 # Access the processed PRS objects
 #head(result$add_prs)  # Processed and scaled additive PRS
 #head(result$int_prs)  # Processed and scaled interaction PRS
 #head(result$cov_prs)  # Processed and scaled covariate PRS
 
 # Check temporary files information
 #if(result$temp_files_saved) {
   #cat("Temporary directory used:", result$temp_dir, "\n")
   #cat("Temporary files created:", paste(result$temp_files, collapse = ", "), "\n")
# }
# Select rows where the third column is NA
na_rows <- add[is.na(add[, 3]), ]
#/data/alh-admzw/plink2 --bfile dummydata --exclude bb.txt --make-bed --out DummyData
