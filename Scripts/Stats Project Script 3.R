# ================================================================
# SCRIPT 3 — Single Index Model
# ================================================================

library(reshape2)
library(ggplot2)
library(readxl)
library(dplyr)

# ------------------------------------------------------------
# 1. Load filtered log returns
# ------------------------------------------------------------

log_returns_df <- readRDS("/Users/hannahguay/Desktop/ECON Final/log_returns_df.rds")
final_tickers  <- read.csv("/Users/hannahguay/Desktop/ECON Final/tickers_filtered_final.csv")$x

log_returns_filtered <- log_returns_df[, c("date", final_tickers)]

# ------------------------------------------------------------
# Define the market factor (equal-weight market)
# ------------------------------------------------------------
market_returns <- rowMeans(log_returns_filtered[, final_tickers], na.rm = TRUE)

# ------------------------------------------------------------
# Run SIM regressions
# ------------------------------------------------------------

alphas <- c()
betas  <- c()
resid_var <- c()

for (tic in final_tickers) {
  
  fit <- lm(log_returns_filtered[[tic]] ~ market_returns)
  
  alphas[tic] <- coef(fit)[1]
  betas[tic]  <- coef(fit)[2]
  resid_var[tic] <- var(residuals(fit))
}

# ------------------------------------------------------------
# Build the SIM covariance matrix
# ------------------------------------------------------------

beta_vec <- matrix(betas, ncol=1)
sigma_m  <- var(market_returns)

sigma_SIM <- beta_vec %*% t(beta_vec) * sigma_m + diag(resid_var)

colnames(sigma_SIM) <- final_tickers
rownames(sigma_SIM) <- final_tickers

write.csv(sigma_SIM, "/Users/hannahguay/Desktop/ECON Final/sigma_SIM.csv")
saveRDS(sigma_SIM, "/Users/hannahguay/Desktop/ECON Final/sigma_SIM.rds")

# Copy full matrix
sigma_SIM_lower <- sigma_SIM

# Replace upper triangle + diagonal with NA
sigma_SIM_lower[upper.tri(sigma_SIM_lower, diag = TRUE)] <- NA

# Export lower triangular CSV
write.csv(
  sigma_SIM_lower,
  "/Users/hannahguay/Desktop/ECON Final/sigma_SIM_lower.csv",
  row.names = TRUE
)

# VISUALS # 

# ------------------------------------------------------------
# Load SIM covariance matrix
# ------------------------------------------------------------
sigma_SIM <- readRDS("/Users/hannahguay/Desktop/ECON Final/sigma_SIM.rds")

sigma_melt <- melt(sigma_SIM, varnames = c("Stock1", "Stock2"), value.name = "Covariance")

# ------------------------------------------------------------
# Heatmap (FULL DATA)
# ------------------------------------------------------------
ggplot(sigma_melt, aes(x = Stock1, y = Stock2, fill = Covariance)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) + 
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    title = "SIM Variance-Covariance Matrix Heatmap",
    x = "",
    y = "",
    fill = "Covariance"
  )

# ------------------------------------------------------------
# Heatmap (LOWER TRIANGLE DATA)
# ------------------------------------------------------------
# Convert lower triangle to long format
sigma_long <- melt(sigma_SIM_lower, varnames = c("Var1", "Var2"), value.name = "value")

# Keep only rows where value is not NA
sigma_long <- sigma_long[!is.na(sigma_long$value), ]

# Plot heatmap
ggplot(sigma_long, aes(x = Var2, y = Var1, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "white") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +   # <-- flips heatmap vertically
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    axis.text.y = element_text(size = 6)
  ) +
  labs(
    title = "SIM Covariance — Lower Triangle",
    x = "",
    y = ""
  )

