# ================================================================
# SCRIPT 5 — Shrinkage VAR Covariance Matrix + Heatmap
# ================================================================

library(dplyr)
library(vars)       # For VAR models
library(corpcor)    # For cov.shrink
library(reshape2)
library(ggplot2)

# ------------------------------------------------------------
# 1. Load log returns + filtered tickers
# ------------------------------------------------------------
log_returns_df <- readRDS("data/processed/log_returns_df.rds")

# Use the same filtered universe as in Script 2
final_tickers <- read.csv("data/processed/tickers_filtered_final.csv",
                          stringsAsFactors = FALSE)[[1]]

# Keep only date + filtered tickers
log_returns_filtered <- log_returns_df %>%
  select(date, all_of(final_tickers))

# Drop any rows with NAs in the selected assets
returns_mat <- log_returns_filtered %>%
  select(-date) %>%
  as.matrix()

complete_idx <- complete.cases(returns_mat)
returns_mat  <- returns_mat[complete_idx, ]

# Create a time series object (monthly frequency assumed)
returns_ts <- ts(returns_mat, frequency = 12)

cat("Number of observations used in VAR:", nrow(returns_mat), "\n")
cat("Number of series (tickers):", ncol(returns_mat), "\n")

# ------------------------------------------------------------
# 2. Select VAR lag order (simple AIC-based choice)
# ------------------------------------------------------------
lag_sel <- VARselect(returns_ts, lag.max = 3, type = "const")
p_opt   <- lag_sel$selection["AIC(n)"]

cat("Selected VAR lag order (AIC):", p_opt, "\n")

# ------------------------------------------------------------
# 3. Estimate VAR model and get residuals
# ------------------------------------------------------------
var_fit <- VAR(returns_ts, p = p_opt, type = "const")

# Matrix of residuals (T x N)
res_mat <- residuals(var_fit)
res_mat <- as.matrix(res_mat)

# ------------------------------------------------------------
# 4. Shrinkage covariance matrix of VAR residuals
#    (Ledoit–Wolf style, via corpcor::cov.shrink)
# ------------------------------------------------------------
sigma_VAR_shrink <- cov.shrink(res_mat)

# Ensure names are aligned with tickers
colnames(sigma_VAR_shrink) <- final_tickers
rownames(sigma_VAR_shrink) <- final_tickers

# ------------------------------------------------------------
# 5. Save shrinkage VAR covariance matrix
# ------------------------------------------------------------
dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/objects", recursive = TRUE, showWarnings = FALSE)

write.csv(
  sigma_VAR_shrink,
  "outputs/tables/sigma_VAR_shrink_full.csv",
  row.names = TRUE
)

saveRDS(
  sigma_VAR_shrink,
  "outputs/objects/sigma_VAR_shrink.rds"
)

cat("Shrinkage VAR covariance matrix saved to:\n",
    "  outputs/tables/sigma_VAR_shrink_full.csv\n")

# ------------------------------------------------------------
# 6. Heatmap of Shrinkage VAR Covariance Matrix
# ------------------------------------------------------------
sigma_VAR_long <- melt(
  sigma_VAR_shrink,
  varnames   = c("Stock1", "Stock2"),
  value.name = "Covariance"
)

ggplot(sigma_VAR_long, aes(x = Stock1, y = Stock2, fill = Covariance)) +
  geom_tile() +
  scale_x_discrete(position = "top") +
  scale_fill_gradient2(
    low      = "blue",
    mid      = "white",
    high     = "red",
    midpoint = 0
  )
