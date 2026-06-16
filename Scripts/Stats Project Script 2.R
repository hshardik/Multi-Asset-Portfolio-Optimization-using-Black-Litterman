# ================================================================
# SCRIPT 2 — Correlation Filtering + Filtered Var/Cov
# ================================================================

library(dplyr)
library(igraph)

# ------------------------------------------------------------
# 1. Load log returns from Script 1
# ------------------------------------------------------------
log_returns_df <- readRDS("data/processed/log_returns_df.rds")

log_returns_mat <- log_returns_df %>%
  select(-date) %>%
  as.matrix()

# ------------------------------------------------------------
# 2. Full Correlation Matrix
# ------------------------------------------------------------
cor_full <- readRDS("outputs/objects/cor_matrix.rds")

# ------------------------------------------------------------
# 3. Correlation Statistics (print to console)
# ------------------------------------------------------------
cor_values <- cor_full[lower.tri(cor_full)]

avg_cor <- mean(cor_values)
max_cor <- max(cor_values)
min_cor <- min(cor_values)

cat("Average Correlation (before filtering): ", avg_cor, "\n")
cat("Largest Correlation (before filtering): ", max_cor, "\n")
cat("Smallest Correlation (before filtering):", min_cor, "\n")

# ------------------------------------------------------------
# 4. Stronger Filtering: remove pairs with |cor| > 0.5
# ------------------------------------------------------------
threshold <- 0.5
adj_matrix <- abs(cor_full) > threshold
diag(adj_matrix) <- FALSE

g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

# Greedy Independent Set Algorithm
greedy_independent_set <- function(graph) {
  Vnames <- V(graph)$name
  remaining <- Vnames
  selected <- c()
  
  while (length(remaining) > 0) {
    v <- remaining[1]
    selected <- c(selected, v)
    neighbors_v <- neighbors(graph, v)
    remaining <- setdiff(remaining, c(v, names(neighbors_v)))
  }
  
  return(selected)
}

filtered_tickers <- greedy_independent_set(g)

cat("Number of securities AFTER filtering:", length(filtered_tickers), "\n")

write.csv(filtered_tickers,
          "data/processed/tickers_filtered.csv",
          row.names = FALSE)

# ------------------------------------------------------------
# 5. Filtered correlation matrix
# ------------------------------------------------------------
log_returns_filtered <- log_returns_mat[, filtered_tickers]

# Remove columns with NA or zero variance
log_returns_filtered <- log_returns_filtered[, colSums(is.na(log_returns_filtered)) == 0]
log_returns_filtered <- log_returns_filtered[, apply(log_returns_filtered, 2, sd) > 0]

cor_filtered <- cor(log_returns_filtered)

write.csv(cor_filtered,
          "outputs/tables/correlation_filtered.csv",
          row.names = TRUE)

saveRDS(cor_filtered, "outputs/objects/cor_filtered.rds")

# New correlation stats
cor_values_filtered <- cor_filtered[lower.tri(cor_filtered)]
avg_cor_after <- mean(cor_values_filtered)

cat("Average Correlation (after filtering): ", avg_cor_after, "\n")

# ------------------------------------------------------------
# 5.2 Update filtered ticker list AFTER cleaning
# ------------------------------------------------------------
final_filtered_tickers <- colnames(log_returns_filtered)

write.csv(final_filtered_tickers,
          "data/processed/tickers_filtered_final.csv",
          row.names = FALSE)

cat("Final number of usable tickers:", length(final_filtered_tickers), "\n")

# ------------------------------------------------------------
# 6. Filtered Variance-Covariance Matrix (+ eigenvalue diagnostics)
# ------------------------------------------------------------
sigma_sample <- cov(log_returns_filtered)

eigen_vals <- eigen(sigma_sample)$values
min_eigen <- min(eigen_vals)
max_eigen <- max(eigen_vals)
condition_number <- max_eigen / min_eigen

cat("\n--- Eigenvalue Diagnostics ---\n")
cat("Min eigenvalue: ", min_eigen, "\n")
cat("Max eigenvalue: ", max_eigen, "\n")
cat("Condition number: ", condition_number, "\n")

# Save as RDS
saveRDS(sigma_sample, "outputs/objects/sigma_sample.rds")

# Export full filtered Var-Cov matrix as CSV
write.csv(sigma_sample,
          "outputs/tables/sigma_sample_full.csv",
          row.names = TRUE)

# ------------------------------------------------------------
# 7. Lower-triangular covariance matrix (for reporting)
# ------------------------------------------------------------
sigma_sample_lower <- sigma_sample
sigma_sample_lower[upper.tri(sigma_sample_lower, diag = TRUE)] <- NA

write.csv(sigma_sample_lower,
          "outputs/tables/sigma_sample_lower.csv",
          row.names = TRUE)

cat("\nScript 2 completed successfully.\n")
