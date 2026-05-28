# Read the file
x <- read.table("Bexp_tar.txt", header = FALSE)

# Rename columns
colnames(x) <- c("FID", "IID", "trait")

# Create binary phenotype with 30% cases
# Top 30% values = cases (1), remaining 70% = controls (0)

threshold <- quantile(x$trait, probs = 0.70, na.rm = TRUE)

x$trait_binary <- ifelse(x$trait > threshold, 1, 0)

# Check prevalence
mean(x$trait_binary)
table(x$trait_binary)

# Keep only binary phenotype
x_binary <- x[, c("FID", "IID", "trait_binary")]

# Save
write.table(
  x_binary,
  "Bexp_tar.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

# Preview
head(x_binary)

































# Read the file
x <- read.table("Bexp_disc_cov.txt", header = FALSE)

# Recode ONLY the 3rd column into binary with 30% cases
# Top 30% values -> 1 (cases)
# Remaining 70% -> 0 (controls)

threshold <- quantile(x$V3, probs = 0.70, na.rm = TRUE)

x$V3 <- ifelse(x$V3 > threshold, 2, 1)

# Check prevalence
mean(x$V3)
table(x$V3)

# Save recoded file
write.table(
  x,
  "Bexp_disc_cov.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)

# Preview
head(x)