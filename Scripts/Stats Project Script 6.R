# ================================================================
# SCRIPT 6 — Efficient Frontier + Optimal Portfolios (vs Frontier)
# ================================================================

library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

# ------------------------------------------------------------
# 0. Helper: safe readRDS or CSV if it exists
# ------------------------------------------------------------
read_matrix_if_exists <- function(rds_path = NULL, csv_path = NULL) {
  if (!is.null(rds_path) && file.exists(rds_path)) {
    mat <- readRDS(rds_path)
    return(as.matrix(mat))
  }
  if (!is.null(csv_path) && file.exists(csv_path)) {
    mat <- read.csv(csv_path, row.names = 1, check.names = FALSE)
    return(as.matrix(mat))
  }
  return(NULL)
}

# ------------------------------------------------------------
# 1. Load log returns + filtered tickers
# ------------------------------------------------------------

# From Script 1
log_returns_df <- readRDS("/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Processed/log_returns_df.rds")

# From Script 2
final_tickers <- read.csv(
  "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/tickers_filtered_final.csv",
  stringsAsFactors = FALSE,
  header = FALSE
)[[1]]

# Some versions of your workflow have a synthetic column "x" (e.g., a portfolio).
# We want only actual asset tickers for the frontier.
asset_tickers <- setdiff(final_tickers, "x")

log_returns_assets <- log_returns_df %>%
  dplyr::select(date, all_of(asset_tickers))

# Matrix of returns (T x N), dropping date
returns_mat <- log_returns_assets %>%
  dplyr::select(-date) %>%
  as.matrix()

# ------------------------------------------------------------
# 2. Mean returns and sample covariance (baseline Σ)
# ------------------------------------------------------------

mu_vec <- colMeans(returns_mat, na.rm = TRUE)  # E[r_i]
names(mu_vec) <- asset_tickers

# Prefer the sigma_sample object from Script 2 (if saved as RDS),
# else fall back to the CSV.
Sigma_sample <- read_matrix_if_exists(
  rds_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/sigma_sample.rds",
  csv_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/sigma_sample_full.csv"
)

if (is.null(Sigma_sample)) {
  stop("Could not find sigma_sample.rds or sigma_sample_full.csv. 
       Run Script 2 before Script 6.")
}

# Align rows/cols with asset_tickers explicitly
Sigma_sample <- Sigma_sample[asset_tickers, asset_tickers]

# Basic dimensions check
N <- length(asset_tickers)
cat("Number of assets in efficient frontier:", N, "\n")

# ------------------------------------------------------------
# 3. Markowitz efficient frontier (no short-sale constraint)
#    Using closed-form formulas (A,B,C,D)
# ------------------------------------------------------------

one_vec   <- rep(1, N)
Sigma_inv <- solve(Sigma_sample)

A <- as.numeric(t(one_vec) %*% Sigma_inv %*% one_vec)
B <- as.numeric(t(one_vec) %*% Sigma_inv %*% mu_vec)
C <- as.numeric(t(mu_vec) %*% Sigma_inv %*% mu_vec)
D <- A * C - B^2

# Global minimum-variance portfolio (GMV) stats
mu_gmv   <- B / A
var_gmv  <- 1 / A
sigma_gmv <- sqrt(var_gmv)

cat("GMV portfolio: E[r] =", mu_gmv, ", sigma =", sigma_gmv, "\n")

# Target-return grid for the efficient frontier
mu_min <- min(mu_vec)
mu_max <- max(mu_vec)

# You can adjust the grid resolution if you want more/less points
target_returns <- seq(mu_min, mu_max, length.out = 60)

frontier_list <- lapply(target_returns, function(mu_target) {
  # Lagrange multipliers for constraints:
  # 1'w = 1, mu'w = mu_target
  lambda1 <- (C - B * mu_target) / D
  lambda2 <- (A * mu_target - B) / D
  
  # Portfolio weights
  w <- Sigma_inv %*% (lambda1 * one_vec + lambda2 * mu_vec)
  w <- as.numeric(w)
  names(w) <- asset_tickers
  
  # Sanity checks
  mu_p <- sum(w * mu_vec)
  var_p <- as.numeric(t(w) %*% Sigma_sample %*% w)
  sigma_p <- sqrt(var_p)
  
  data.frame(
    target_return = mu_target,
    realized_return = mu_p,
    variance = var_p,
    sigma = sigma_p
  )
})

frontier_df <- bind_rows(frontier_list)

# ------------------------------------------------------------
# 4. Export efficient frontier table
# ------------------------------------------------------------

#dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)
#dir.create("outputs/objects", recursive = TRUE, showWarnings = FALSE)
#dir.create("outputs/plots",  recursive = TRUE, showWarnings = FALSE)

write.csv(
  frontier_df,
  "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/efficient_frontier.csv",
  row.names = FALSE
)

saveRDS(
  frontier_df,
  "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/efficient_frontier.rds"
)

cat("Efficient frontier table saved to outputs/tables/efficient_frontier.csv\n")

# ------------------------------------------------------------
# 5. Asset-level risk/return table (points on the graph)
# ------------------------------------------------------------

asset_stats <- data.frame(
  ticker      = asset_tickers,
  exp_return  = as.numeric(mu_vec),
  sigma       = sqrt(diag(Sigma_sample))
)

write.csv(
  asset_stats,
  "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/asset_risk_return.csv",
  row.names = FALSE
)

cat("Asset risk/return table saved to outputs/tables/asset_risk_return.csv\n")

# ------------------------------------------------------------
# 6. “Optimal” portfolios using different covariance estimators
#    (Sample, Constant Correlation, SIM, Shrinkage), if available
# ------------------------------------------------------------

risk_free <- 0.0003  # monthly risk-free rate used in Script 1; adjust if needed

tangency_portfolio <- function(mu, Sigma, rf) {
  excess <- mu - rf
  Sigma_inv <- solve(Sigma)
  w_raw <- Sigma_inv %*% excess
  w <- as.numeric(w_raw / sum(w_raw))
  names(w) <- names(mu)
  
  mu_p <- sum(w * mu)
  var_p <- as.numeric(t(w) %*% Sigma %*% w)
  sigma_p <- sqrt(var_p)
  
  list(weights = w, mu = mu_p, sigma = sigma_p)
}

cov_list <- list(
  Sample = Sigma_sample
)

# Constant-correlation covariance (Script 4)
Sigma_const <- read_matrix_if_exists(
  rds_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/constant_correlation_covariance_matrix.rds",
  csv_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/constant_correlation_covariance_matrix.csv"
)
if (!is.null(Sigma_const)) {
  Sigma_const <- Sigma_const[asset_tickers, asset_tickers]
  cov_list[["ConstantCorr"]] <- Sigma_const
}

# SIM covariance (Script 3)
Sigma_SIM <- read_matrix_if_exists(
  rds_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/sigma_SIM.rds",
  csv_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/sigma_SIM.csv"   # adjust if your Script 3 saves elsewhere
)
if (!is.null(Sigma_SIM)) {
  Sigma_SIM <- Sigma_SIM[asset_tickers, asset_tickers]
  cov_list[["SIM"]] <- Sigma_SIM
}

# Shrinkage covariance (Script 5)
Sigma_shrink <- read_matrix_if_exists(
  rds_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/sigma_shrinkage.rds",
  csv_path = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/sigma_shrinkage_full.csv"
)
if (!is.null(Sigma_shrink)) {
  Sigma_shrink <- Sigma_shrink[asset_tickers, asset_tickers]
  cov_list[["Shrinkage"]] <- Sigma_shrink
}

# Compute tangency portfolios for each estimator
opt_list <- list()
weights_long_list <- list()

for (nm in names(cov_list)) {
  Sigma_est <- cov_list[[nm]]
  tp <- tangency_portfolio(mu_vec, Sigma_est, risk_free)
  
  opt_list[[nm]] <- data.frame(
    method      = nm,
    exp_return  = tp$mu,
    sigma       = tp$sigma
  )
  
  weights_long_list[[nm]] <- data.frame(
    method = nm,
    ticker = names(tp$weights),
    weight = tp$weights
  )
}

optimal_ports_df <- bind_rows(opt_list)
optimal_weights_df <- bind_rows(weights_long_list)

write.csv(
  optimal_ports_df,
  "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/optimal_portfolios_summary.csv",
  row.names = FALSE
)

write.csv(
  optimal_weights_df,
  "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Tables/optimal_portfolios_weights.csv",
  row.names = FALSE
)

cat("Optimal portfolio summaries saved to outputs/tables/optimal_portfolios_summary.csv\n")
cat("Optimal portfolio weights saved to outputs/tables/optimal_portfolios_weights.csv\n")

# ------------------------------------------------------------
# 7. Plot: Efficient Frontier + Individual Assets + Optimal Portfolios
# ------------------------------------------------------------

# Frontier: risk (sigma) on x-axis, expected return on y-axis
p <- ggplot() +
  # Efficient frontier (line)
  geom_line(
    data = frontier_df,
    aes(x = sigma, y = realized_return),
    linewidth = 1,
    color = "steelblue"
  ) +
  # Individual assets (red dots)
  geom_point(
    data = asset_stats,
    aes(x = sigma, y = exp_return),
    color = "red",
    size = 2,
    alpha = 0.7
  ) +
  # Optimal portfolios from each covariance estimator (triangles)
  geom_point(
    data = optimal_ports_df,
    aes(x = sigma, y = exp_return, color = method),
    size = 3,
    shape = 17
  ) +
  labs(
    title = "Efficient Frontier and Optimal Portfolios",
    x = "Portfolio Volatility (σ)",
    y = "Expected Return",
    color = "Estimator"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave(
  filename = "/Users/hardik/Documents/Business School/Intro to Prob and Stats/Project/Outputs/Plots/efficient_frontier_with_optimal_portfolios.png",
  plot     = p,
  width    = 8,
  height   = 5,
  dpi      = 300
)

cat("Plot saved to outputs/plots/efficient_frontier_with_optimal_portfolios.png\n")
cat("Script 6 completed successfully.\n")

