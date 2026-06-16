# ================================================================
# SCRIPT 4 — Constant Correlation Matrix (+ optional covariance)
# ================================================================

library(dplyr)
library(reshape2)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load filtered correlation matrix (for tickers & dimension)
# ------------------------------------------------------------
# Path assumes Script 2's output structure; adjust if needed.
cor_filtered <- read.csv(
  "outputs/tables/correlation_filtered.csv",
  row.names   = 1,
  check.names = FALSE
)

tickers <- colnames(cor_filtered)
n       <- length(tickers)

cat("Number of filtered securities:", n, "\n")

# ------------------------------------------------------------
# 2. Set the average (constant) correlation
# ------------------------------------------------------------
avg_cor <- 0.2125   # given

# ------------------------------------------------------------
# 3. Build the Constant Correlation Matrix
#    - diagonal = 1
#    - off-diagonal = avg_cor
# ------------------------------------------------------------
const_cor_mat <- matrix(avg_cor, nrow = n, ncol = n)
diag(const_cor_mat) <- 1

rownames(const_cor_mat) <- tickers
colnames(const_cor_mat) <- tickers

# ------------------------------------------------------------
# 4. Export Constant Correlation Matrix as CSV
# ------------------------------------------------------------
write.csv(
  const_cor_mat,
  "outputs/tables/constant_correlation_matrix.csv",
  row.names = TRUE
)

cat("Constant correlation matrix saved to outputs/tables/constant_correlation_matrix.csv\n")

# ------------------------------------------------------------
# 4b. OPTIONAL: Build Constant-Correlation Covariance Matrix
#     (analogous to your VBA function using variances)
#     Sigma_const = D * C_const * D, where D = diag(sigma_i)
# ------------------------------------------------------------
# Only run this if you also want the covariance estimator
if (file.exists("outputs/tables/sigma_sample_full.csv")) {
  
  sigma_sample <- read.csv(
    "outputs/tables/sigma_sample_full.csv",
    row.names   = 1,
    check.names = FALSE
  )
  
  # Align ordering with 'tickers' to be safe
  sigma_sample <- as.matrix(sigma_sample[tickers, tickers])
  
  # Standard deviations from sample covariance
  sd_vec <- sqrt(diag(sigma_sample))
  D      <- diag(sd_vec)
  
  # Constant-correlation covariance matrix
  const_cov_mat <- D %*% const_cor_mat %*% D
  
  rownames(const_cov_mat) <- tickers
  colnames(const_cov_mat) <- tickers
  
  write.csv(
    const_cov_mat,
    "outputs/tables/constant_correlation_covariance_matrix.csv",
    row.names = TRUE
  )
  
  cat("Constant-correlation covariance matrix saved to outputs/tables/constant_correlation_covariance_matrix.csv\n")
}

# ------------------------------------------------------------
# 5. Heatmap of Constant Correlation Matrix
# ------------------------------------------------------------
const_cor_long <- melt(
  const_cor_mat,
  varnames  = c("Stock1", "Stock2"),
  value.name = "Correlation"
)

ggplot(const_cor_long, aes(x = Stock1, y = Stock2, fill = Correlation)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  scale_fill_gradient2(
    low = "blue", mid = "white", high = "red",
    midpoint = 0
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    title = "Constant Correlation Matrix (ρ̄ = 0.2125)",
    x = "",
    y = "",
    fill = "Correlation"
  )

cat("Script 4 completed successfully.\n")
