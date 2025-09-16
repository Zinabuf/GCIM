## For Continuous Outcomes
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "/data/alh-admzw/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qphen_disc.txt", "Qcov_disc.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qphen_disc_cov.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe) 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)

# # Step 4: Run GCIM analysis with automatic saving and scaling
 result0 <- gcim_q0("Qexp_tar.txt", "Qphen_tar_cov.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 result1 <- gcim_q1("Qexp_tar.txt", "Qphen_tar_cov.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)

 
 # Step 5: Access results and processed PRS objects
 print(result0$model_summary)

 print(result1$model_summary)

 cat("Sample size:", result$sample_size, "\n")
cat("R-squared:", result0$r_squared, "\n")
 cat("R-squared:", result1$r_squared, "\n")
###############################################
cat("Adjusted R-squared:", result0$adj_r_squared, "\n")
 cat("Adjusted R-squared:", result1$adj_r_squared, "\n")

 cat("Residual standard error:", result$residual_se, "\n")
 cat("F-statistic:", result0$f_statistic[1], "on", result0$f_statistic[2], "and", result0$f_statistic[3], "DF\n")
 cat("AIC:", result0$aic, "\n")
 cat("BIC:", result0$bic, "\n")
 
 # Access the processed PRS objects
 #head(result$add_prs)  # Processed and scaled additive PRS
 #head(result$int_prs)  # Processed and scaled interaction PRS
 #head(result$cov_prs)  # Processed and scaled covariate PRS
 
 # Check temporary files information
 #if(result$temp_files_saved) {
  # cat("Temporary directory used:", result$temp_dir, "\n")
  # cat("Temporary files created:", paste(result$temp_files, collapse = ", "), "\n")
 #}
 
 # Additional model diagnostics
 # Check residuals
 plot(result$model)  # Diagnostic plots
 
 # Extract coefficients with confidence intervals
 coef_table <- stats::confint(result$model)
 print(coef_table)