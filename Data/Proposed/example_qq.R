## For Continuous Outcomes
# Load required libraries
library(GxEprs)
library(GCIM)
# Set plink path
plink_path <- "/data/alh-admzw/plink2"

# For quantitative traits, use corresponding functions
a <- GWAS_quantitative(plink_path, "DummyData", "Qexp_disc.txt", "Qcov_disc.txt")
b <- GWEIS_quantitative(plink_path, "DummyData", "Qphen_disc.txt", "Qexp_dis_cov.txt")

# Extract and compute PRS
trd <- a[c("ID", "A1", "BETA")]
add <- b[c("ID", "A1", "ADD_BETA")]
gxe <- b[c("ID", "A1", "INTERACTION_BETA")]

q <- PRS_quantitative(plink_path, "DummyData", summary_input = add)
write.table(q, "add_prs.txt", row.names=F, quote = F, sep = "\t")
r <- PRS_quantitative(plink_path, "DummyData", summary_input = gxe)
write.table(r, "int_prs.txt", row.names=F, quote = F, sep = "\t") 
p <- PRS_quantitative(plink_path, "DummyData", summary_input = trd)
write.table(p, "cov_prs.txt", row.names=F, quote = F, sep = "\t")


# # Step 4: Run GCIM analysis with automatic saving and scaling
 result1 <- gcim_q0("Qphen_tar.txt", "Qexp_tar_cov.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)
 result2 <- gcim_q1("Qphen_tar.txt", "Qexp_tar_cov.txt", 
                  Add_PRS = q, Int_PRS = r, Cov_PRS = p)

 # Step 5: Access results and processed PRS objects
 print(result1$model_summary)
 print(result2$model_summary)
###############################################
 cat("Sample size:", result1$sample_size, "\n")
 cat("R-squared:", result1$r_squared, "\n")
 cat("Adjusted R-squared:", result1$adj_r_squared, "\n")
 cat("Residual standard error:", result1$residual_se, "\n")
 cat("F-statistic:", result1$f_statistic[1], "on", result1$f_statistic[2], "and", result1$f_statistic[3], "DF\n")
 cat("AIC:", result1$aic, "\n")
 cat("BIC:", result1$bic, "\n")
 
