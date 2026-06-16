# ================================================================
# SCRIPT 7 — Black-Litterman Benchmark + Optimal Portfolio
# ================================================================

library(dplyr)
library(ggplot2)
library(reshape2)

# ------------------------------------------------------------
# 0. Helper: read matrix from RDS or CSV if it exists
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

dir.create("outputs/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/objects", recursive = TRUE, showWarnings = FALSE)
dir.create("outputs/plots",   recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 1. Load log returns + filtered ticker universe
# ------------------------------------------------------------

log_returns_df <- readRDS("data/processed/log_returns_df.rds")

final_tickers <- read.csv(
  "data/processed/tickers_filtered_final.csv",
  stringsAsFactors = FALSE,
  header = FALSE
)[[1]]

asset_tickers <- final_tickers
N <- length(asset_tickers)

log_returns_assets <- log_returns_df %>%
  dplyr::select(date, all_of(asset_tickers))

returns_mat <- as.matrix(log_returns_assets[, -1])

# Sample mean returns
mu_sample <- colMeans(returns_mat, na.rm = TRUE)
names(mu_sample) <- asset_tickers

# ------------------------------------------------------------
# 2. Choose covariance matrix Σ
#    (prefer shrinkage; fall back to sample Σ)
# ------------------------------------------------------------

Sigma <- read_matrix_if_exists(
  rds_path = "outputs/objects/sigma_shrinkage.rds",
  csv_path = "outputs/tables/sigma_shrinkage_full.csv"
)

if (is.null(Sigma)) {
  Sigma <- read_matrix_if_exists(
    rds_path = "outputs/objects/sigma_sample.rds",
    csv_path = "outputs/tables/sigma_sample_full.csv"
  )
}

if (is.null(Sigma)) {
  Sigma <- cov(returns_mat, use = "pairwise.complete.obs")
}

Sigma <- Sigma[asset_tickers, asset_tickers]
one_vec <- rep(1, N)

# ------------------------------------------------------------
# 3. Benchmark portfolio weights (like “Benchmark proportions”)
#    - If data/inputs/benchmark_weights.csv exists, use it
#      (columns: ticker, weight).
#    - Otherwise, default to equal-weight benchmark.
# ------------------------------------------------------------

bench_path <- "data/inputs/benchmark_weights.csv"

if (file.exists(bench_path)) {
  bench_df <- read.csv(bench_path, stringsAsFactors = FALSE)
  w_bench <- rep(0, N); names(w_bench) <- asset_tickers
  w_bench[bench_df$ticker] <- bench_df$weight
} else {
  w_bench <- rep(1 / N, N)
  names(w_bench) <- asset_tickers
}

# Normalize just in case
w_bench <- w_bench / sum(w_bench)

# ------------------------------------------------------------
# 4. Benchmark expected return, variance, and normalizing factor δ
#    (mimics “Anticipated benchmark return”, “Current T-bill rate”,
#     and “Normalizing factor” in the Excel sheet)
# ------------------------------------------------------------

risk_free <- 0.0003   # monthly risk-free rate (adjust if desired)

bench_returns <- as.numeric(returns_mat %*% w_bench)
mu_bench <- mean(bench_returns, na.rm = TRUE)

var_bench <- as.numeric(t(w_bench) %*% Sigma %*% w_bench)
sigma_bench <- sqrt(var_bench)

# Normalizing factor:
#   δ = (E[R_bench] - r_f) / Var(R_bench)
delta <- (mu_bench - risk_free) / var_bench

cat("Benchmark expected return (data-based):", mu_bench, "\n")
cat("Benchmark volatility:", sigma_bench, "\n")
cat("Normalizing factor (delta):", delta, "\n")

# ------------------------------------------------------------
# 5. Equilibrium (implied) returns — prior for Black-Litterman
#    π = δ Σ w_bench   (excess returns)
#    μ_eq = r_f + π    (total returns)
#    This is the “with normalizing factor” column in Excel.
# ------------------------------------------------------------

pi_excess <- as.numeric(delta * (Sigma %*% w_bench))  # excess
names(pi_excess) <- asset_tickers

mu_equil <- risk_free + pi_excess

equil_df <- data.frame(
  ticker           = asset_tickers,
  implied_excess   = pi_excess,
  implied_total    = mu_equil
)

write.csv(
  equil_df,
  "outputs/tables/bl_equilibrium_returns.csv",
  row.names = FALSE
)
saveRDS(
  equil_df,
  "outputs/objects/bl_equilibrium_returns.rds"
)

# ------------------------------------------------------------
# 6. Define Black-Litterman views (P, Q, Ω)
#    Here we create ONE illustrative view:
#    - The asset with the highest sample mean return is expected
#      to earn its sample mean (excess) going forward.
#    You can edit this block to match your own views.
# ------------------------------------------------------------

# Asset with highest sample mean
view_asset <- names(sort(mu_sample, decreasing = TRUE))[1]
cat("View asset (highest sample mean):", view_asset, "\n")

K <- 1  # number of views

P <- matrix(0, nrow = K, ncol = N,
            dimnames = list(paste0("View", 1:K), asset_tickers))
P[1, which(asset_tickers == view_asset)] <- 1

# Q is the view on EXCESS return of that asset
Q <- matrix(mu_sample[view_asset] - risk_free,
            nrow = K, ncol = 1,
            dimnames = list(paste0("View", 1:K), "Q"))

# Tau: scaling of prior uncertainty
tau <- 0.05   # common choice; you can adjust

# Omega: view uncertainty matrix (diagonal)
# Here we use: Ω = diag( diag(P (τΣ) P') )
Omega <- diag(diag(P %*% (tau * Sigma) %*% t(P)))
rownames(Omega) <- colnames(Omega) <- rownames(P)

# Export views
write.csv(P, "outputs/tables/bl_views_P.csv")
write.csv(Q, "outputs/tables/bl_views_Q.csv")
write.csv(Omega, "outputs/tables/bl_views_Omega.csv")

# ------------------------------------------------------------
# 7. Black-Litterman posterior (μ_BL, Σ_BL)
#    μ_BL (excess returns) =
#      [ (τΣ)^(-1) + P' Ω^(-1) P ]^(-1) [ (τΣ)^(-1) π + P' Ω^(-1) Q ]
#    Σ_BL = Σ + [ (τΣ)^(-1) + P' Ω^(-1) P ]^(-1)
# ------------------------------------------------------------

tauSigma_inv <- solve(tau * Sigma)
Omega_inv    <- solve(Omega)

M <- tauSigma_inv + t(P) %*% Omega_inv %*% P

mu_bl_excess <- as.numeric(
  solve(M, tauSigma_inv %*% pi_excess + t(P) %*% Omega_inv %*% Q)
)
names(mu_bl_excess) <- asset_tickers

mu_bl_total <- risk_free + mu_bl_excess

Sigma_bl <- Sigma + solve(M)

# Export posterior returns and covariance
bl_returns_df <- data.frame(
  ticker          = asset_tickers,
  bl_excess       = mu_bl_excess,
  bl_total        = mu_bl_total,
  sample_mean     = mu_sample
)

write.csv(
  bl_returns_df,
  "outputs/tables/bl_posterior_returns.csv",
  row.names = FALSE
)
saveRDS(
  Sigma_bl,
  "outputs/objects/bl_posterior_covariance.rds"
)
write.csv(
  Sigma_bl,
  "outputs/tables/bl_posterior_covariance_full.csv",
  row.names = TRUE
)

# ------------------------------------------------------------
# 8. BL optimal portfolio weights
#    (tangency-style: w ∝ Σ_BL^(-1) μ_BL_excess, normalized to sum 1)
# ------------------------------------------------------------

Sigma_bl_inv <- solve(Sigma_bl)

w_bl_unscaled <- Sigma_bl_inv %*% mu_bl_excess
w_bl <- as.numeric(w_bl_unscaled / sum(w_bl_unscaled))
names(w_bl) <- asset_tickers

# Benchmark metrics (using BL expected returns)
mu_bench_bl   <- sum(w_bench * mu_bl_total)
sigma_bench_bl <- sqrt(as.numeric(t(w_bench) %*% Sigma_bl %*% w_bench))

# BL optimal portfolio metrics
mu_bl_port   <- sum(w_bl * mu_bl_total)
sigma_bl_port <- sqrt(as.numeric(t(w_bl) %*% Sigma_bl %*% w_bl))

bl_weights_df <- data.frame(
  ticker = asset_tickers,
  benchmark_weight = w_bench,
  bl_opt_weight    = w_bl
)

write.csv(
  bl_weights_df,
  "outputs/tables/bl_weights_benchmark_and_optimal.csv",
  row.names = FALSE
)

bl_summary_df <- data.frame(
  portfolio = c("Benchmark", "BL_Optimal"),
  exp_return = c(mu_bench_bl, mu_bl_port),
  sigma      = c(sigma_bench_bl, sigma_bl_port)
)

write.csv(
  bl_summary_df,
  "outputs/tables/bl_portfolio_summary.csv",
  row.names = FALSE
)

cat("BL optimal portfolio expected return:", mu_bl_port, "\n")
cat("BL optimal portfolio volatility:", sigma_bl_port, "\n")

# ------------------------------------------------------------
# 9. Efficient frontier (using BL posterior μ & Σ)
# ------------------------------------------------------------

Sigma_inv_bl <- Sigma_bl_inv
A <- as.numeric(t(one_vec) %*% Sigma_inv_bl %*% one_vec)
B <- as.numeric(t(one_vec) %*% Sigma_inv_bl %*% mu_bl_total)
C <- as.numeric(t(mu_bl_total) %*% Sigma_inv_bl %*% mu_bl_total)
D <- A * C - B^2

target_returns <- seq(min(mu_bl_total), max(mu_bl_total), length.out = 60)

frontier_list <- lapply(target_returns, function(mu_t) {
  lambda1 <- (C - B * mu_t) / D
  lambda2 <- (A * mu_t - B) / D
  
  w <- Sigma_inv_bl %*% (lambda1 * one_vec + lambda2 * mu_bl_total)
  w <- as.numeric(w)
  
  mu_p <- sum(w * mu_bl_total)
  var_p <- as.numeric(t(w) %*% Sigma_bl %*% w)
  sigma_p <- sqrt(var_p)
  
  data.frame(
    target_return = mu_t,
    exp_return    = mu_p,
    sigma         = sigma_p
  )
})

frontier_bl_df <- bind_rows(frontier_list)

write.csv(
  frontier_bl_df,
  "outputs/tables/bl_efficient_frontier.csv",
  row.names = FALSE
)

# Asset-level BL risk/return
asset_stats_bl <- data.frame(
  ticker     = asset_tickers,
  exp_return = mu_bl_total,
  sigma      = sqrt(diag(Sigma_bl))
)

write.csv(
  asset_stats_bl,
  "outputs/tables/bl_asset_risk_return.csv",
  row.names = FALSE
)

# ------------------------------------------------------------
# 10. Plot: BL Efficient Frontier + Assets + Benchmark + BL Optimal
# ------------------------------------------------------------

bench_point <- data.frame(
  portfolio  = "Benchmark",
  exp_return = mu_bench_bl,
  sigma      = sigma_bench_bl
)

bl_point <- data.frame(
  portfolio  = "BL_Optimal",
  exp_return = mu_bl_port,
  sigma      = sigma_bl_port
)

p_bl <- ggplot() +
  # BL efficient frontier
  geom_line(
    data = frontier_bl_df,
    aes(x = sigma, y = exp_return),
    linewidth = 1,
    color = "steelblue"
  ) +
  # Individual assets
  geom_point(
    data = asset_stats_bl,
    aes(x = sigma, y = exp_return),
    color = "red",
    size = 2,
    alpha = 0.7
  ) +
  # Benchmark
  geom_point(
    data = bench_point,
    aes(x = sigma, y = exp_return),
    shape = 17,
    size = 3,
    color = "black"
  ) +
  # BL optimal portfolio
  geom_point(
    data = bl_point,
    aes(x = sigma, y = exp_return),
    shape = 17,
    size = 3,
    color = "darkgreen"
  ) +
  labs(
    title = "Black-Litterman Efficient Frontier and Optimal Portfolio",
    x = "Volatility (σ)",
    y = "Expected Return",
    color = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "bottom"
  )

ggsave(
  filename = "outputs/plots/bl_efficient_frontier_with_optimal.png",
  plot     = p_bl,
  width    = 8,
  height   = 5,
  dpi      = 300
)

cat("Script 6 (Black-Litterman) completed successfully.\n")
