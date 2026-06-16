# ================================================================
# SCRIPT 1 — Basic Returns + Var/Cov + Correlation
# ================================================================

library(readxl)
library(dplyr)
library(lubridate)
library(xts)
library(PerformanceAnalytics)

# ------------------------------------------------------------
# 1. Load monthly adjusted price data
# ------------------------------------------------------------
prices_df <- read_excel("project_data.xlsx") %>%
  mutate(date = as_date(date)) %>%
  arrange(date)

prices_xts <- xts(prices_df[,-1], order.by = prices_df$date)

# ------------------------------------------------------------
# 2. Monthly log returns
# ------------------------------------------------------------
log_returns_xts <- Return.calculate(prices_xts, method = "log")
log_returns_xts <- log_returns_xts[-1, ]   # remove NA row

log_returns_df <- data.frame(
  date = index(log_returns_xts),
  coredata(log_returns_xts)
)

# Save log returns for later use
saveRDS(log_returns_df, "data/processed/log_returns_df.rds")

# Matrix (no date)
log_returns_mat <- log_returns_df %>%
  select(-date) %>%
  as.matrix()

# ------------------------------------------------------------
# 3. Excess returns
# ------------------------------------------------------------
risk_free <- 0.0003

excess_returns_df <- log_returns_df %>%
  mutate(across(-date, ~ .x - risk_free))

saveRDS(excess_returns_df, "data/processed/excess_returns_df.rds")

# ------------------------------------------------------------
# 4. Variance-Covariance Matrix
# ------------------------------------------------------------
var_cov_matrix <- cov(log_returns_mat, use = "pairwise.complete.obs")

# Save as RDS
saveRDS(var_cov_matrix, "outputs/objects/var_cov_matrix.rds")

# 🔥 NEW: Export full Var-Cov matrix as CSV
write.csv(var_cov_matrix,
          "outputs/tables/var_cov_matrix_full.csv",
          row.names = TRUE)

# ------------------------------------------------------------
# 5. Correlation Matrix + CSV Export
# ------------------------------------------------------------
cor_matrix <- cor(log_returns_mat, use = "pairwise.complete.obs")

write.csv(cor_matrix,
          "outputs/tables/correlation_full.csv",
          row.names = TRUE)

saveRDS(cor_matrix, "outputs/objects/cor_matrix.rds")

cat("Script 1 completed successfully.\n")
