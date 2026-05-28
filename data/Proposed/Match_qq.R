x <- read.table("Bphen_tar.txt", header = FALSE)
y <- read.table("Bexp_tar_cov.txt", header = FALSE)

# Remove accidental header row from x if present
x <- x[x$V1 != "V1" & x$V2 != "V2", ]

# Convert phenotype to numeric
x$V3 <- as.numeric(x$V3)

# Merge by FID and IID
merged <- merge(y, x[, c("V1", "V2", "V3")],
                by = c("V1", "V2"),
                all.x = TRUE)

# Replace y$V3 with x$V3
merged$V3.x <- merged$V3.y

# Replace y$V4 with (x$V3)^2
#merged$V4 <- merged$V3.y^2

# Remove extra phenotype column
merged$V3.y <- NULL

# Restore column names
colnames(merged)[1:4] <- c("FID", "IID", "V3", "V4")

# Save
write.table(
  merged,
  "Bphen_tar_cov.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = "\t"
)