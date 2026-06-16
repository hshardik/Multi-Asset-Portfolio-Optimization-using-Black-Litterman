# Multi-Asset Portfolio Optimization using Black–Litterman

A full end-to-end portfolio construction pipeline in R — from a raw universe of 138 global equities, ETFs, and funds to a Bayesian-enhanced Black–Litterman optimal portfolio. The project implements and compares multiple covariance estimators, constructs Markowitz efficient frontiers, and derives a benchmark-anchored BL portfolio that blends equilibrium returns with explicit investor views.

---

## Project Objective

The central question this project answers is:

> *How do we go from a big, noisy, overlapping set of names to a stable, diversified optimal portfolio using modern portfolio theory and Bayesian enhancements (Black–Litterman)?*

Starting from a 5-year monthly dataset of 138 securities spanning technology, financials, consumer discretionary, industrials, healthcare, and utilities, the project:

1. Engineers log, excess, and mean returns as inputs to all downstream models.
2. Implements and compares four covariance estimators — sample, Single Index Model (SIM), constant-correlation, and shrinkage — to assess stability and portfolio implications.
3. Solves Markowitz mean–variance optimization to construct efficient frontiers and tangency portfolios under each risk model.
4. Develops a Black–Litterman framework using a consensus benchmark and investor views to derive posterior returns, generate a BL efficient frontier, and compute an optimal BL portfolio.

---

## Security Universe

138 securities across sectors and instrument types (individual stocks, ETFs, sector/thematic funds):

| Sector | Sub-category | Tickers |
|---|---|---|
| Technology | Semiconductors & Equipment | AMD, AVGO, QCOM, TXN, AMAT, KLAC, LRCX, MU, ADI, MCHP, MPWR, ON, SWKS |
| Technology | Software & Cloud | MSFT, CRM, NOW, ADBE, INTU, CRWD, PANW, SNPS, CDNS, WDAY, DDOG, ORCL, FTNT |
| Technology | IT Services & Hardware | IBM, CSCO, AKAM, NTAP, AAPL, DELL, HPE, GLW, FLEX, WDC, STX, SMCI |
| Financials | Investment Banks & Asset Management | GS, MS, BLK, BX, KKR, APO, SCHW, TROW, BK, AMP |
| Financials | Exchanges & Market Infrastructure | CME, ICE, CBOE |
| Financials | Payments & Processing | MA, V, PYPL, FIS, FI |
| Financials | Consumer Finance & Services | COF, SYF, IBKR, VRSK, MSCI |
| Consumer Discretionary | E-Commerce & Retail | AMZN, HD, LOW, TJX, ROST, ULTA, AZO, ORLY, BBY |
| Consumer Discretionary | Restaurants | MCD, CMG, SBUX, YUM, DPZ |
| Consumer Discretionary | Apparel & Footwear | NKE, DECK, RL, CROX |
| Consumer Discretionary | Hotels & Gaming | MAR, HLT, MGM, LVS |
| Consumer Discretionary | Automotive | TSLA, GM |
| Industrials | Aerospace & Defense | GE, TDG |
| Industrials | Construction & Building | NVR, DHI, LEN |
| Industrials | Diversified Industrials & Equipment | ROP, TDY, CPRT, HUBB, MSI, POOL |
| Communication Services | — | EBAY, LYV |
| Healthcare | Pharma & Biotech | LLY, ABT, TMO, AMGN, PHAR, NRIX, ELV, CHGCY, TAK |
| Healthcare | Funds & ETFs | ETIHX, FBIOX, FBTDX, VHCIX, PRHSX, SHSKX, PJP, IXJ, XBI, FXH |
| Utilities | — | NEE, DUK, D, AEP, SO, ED, EIX, EXC, XEL, AES, CNP, ETR, EVRG, LNT, ATO, NI, NRG, AWK, WTRG, PEG, PPL |

---

## Repository Structure

```
.
├── Scripts/
│   ├── Stats Project Script 1.R   # Data preparation & baseline risk metrics
│   ├── Stats Project Script 2.R   # Correlation filtering & universe reduction
│   ├── Stats Project Script 3.R   # Single Index Model (SIM) covariance
│   ├── Stats Project Script 4.R   # Constant correlation covariance
│   ├── Stats Project Script 5.R   # Shrinkage covariance
│   ├── Stats Project Script 6.R   # Markowitz efficient frontier & tangency portfolios
│   └── Stats Project Script 7.R   # Black–Litterman framework & BL optimal portfolio
│
├── Inputs/
│   ├── project_data.xlsx           # Monthly adjusted prices for all 138 tickers (5 years)
│   ├── R Project Stock Price Data.xlsx  # Raw price data workbook
│   └── benchmark_weights.csv       # Consensus benchmark weights (avg of 4 tangency portfolios)
│
├── Outputs/
│   ├── Plots/
│   │   ├── efficient_frontier_with_optimal_portfolios.png  # Markowitz frontier + 4 tangency portfolios
│   │   ├── bl_efficient_frontier_with_optimal.png          # BL frontier + benchmark + BL optimal
│   │   ├── Constant-Corr Matrix.png                        # Heatmap: constant-correlation Σ
│   │   ├── SIM VarCovar Matrix.png                         # Heatmap: SIM covariance matrix
│   │   ├── Single Index Model VarCovar.png                 # Heatmap: SIM covariance (alt view)
│   │   └── Shrinkage-var-Covar.png                         # Heatmap: shrinkage covariance matrix
│   │
│   └── Tables/
│       ├── asset_risk_return.csv               # Per-asset expected return and volatility
│       ├── correlation_full.csv                # Full 138×138 correlation matrix
│       ├── correlation_filtered.csv            # Filtered correlation matrix (post graph reduction)
│       ├── sigma_sample_full.csv               # Sample covariance (filtered universe, full)
│       ├── sigma_sample_lower.csv              # Sample covariance (lower triangle)
│       ├── sigma_SIM.csv                       # SIM covariance matrix
│       ├── sigma_SIM_lower.csv                 # SIM covariance (lower triangle)
│       ├── constant_correlation_matrix.csv     # Constant correlation matrix
│       ├── constant_correlation_covariance_matrix.csv  # Constant-corr covariance Σ
│       ├── sigma_shrinkage_full.csv            # Shrinkage covariance matrix
│       ├── sigma_VAR_shrink_full.csv           # Variance-adjusted shrinkage covariance
│       ├── efficient_frontier.csv              # Markowitz frontier (σ, μ) grid
│       ├── optimal_portfolios_weights.csv      # Tangency weights under all 4 Σ estimators
│       ├── optimal_portfolios_summary.csv      # Return, risk, Sharpe for each tangency portfolio
│       ├── tickers_filtered.csv                # Tickers after correlation filter
│       ├── tickers_filtered_final.csv          # Final cleaned ticker list
│       ├── bl_equilibrium_returns.csv          # BL implied equilibrium returns (π)
│       ├── bl_posterior_returns.csv            # BL posterior excess returns (μ_BL)
│       ├── bl_posterior_covariance_full.csv    # BL posterior covariance (Σ_BL)
│       ├── bl_asset_risk_return.csv            # Per-asset risk/return under BL
│       ├── bl_efficient_frontier.csv           # BL efficient frontier grid
│       ├── bl_portfolio_summary.csv            # BL optimal portfolio return, risk, Sharpe
│       ├── bl_weights_benchmark_and_optimal.csv # Benchmark vs. BL optimal weights
│       ├── bl_views_P.csv                      # View pick matrix (K×N)
│       ├── bl_views_Q.csv                      # View return vector (K×1)
│       └── bl_views_Omega.csv                  # View uncertainty covariance (K×K)
│
├── Project Presentation.pdf        # Slide deck summarizing methodology and results
└── Readme.pdf                      # Detailed technical writeup
```

> **Note:** `Outputs/Objects/` (`.rds` files) are excluded from version control via `.gitignore` as they are regenerable by running the scripts in order.

---

## Methodology

### Script 1 — Data Preparation & Baseline Risk Metrics

Loads monthly adjusted prices from `project_data.xlsx` and converts them to an `xts` time-series object. Computes log returns for all 138 assets:

```
r_{i,t} = ln(P_{i,t} / P_{i,t-1})
```

Constructs excess returns over a monthly risk-free rate (r_f = 0.0003), and estimates the full 138×138 sample variance–covariance matrix Σ = cov(R) and correlation matrix ρ = cor(R).

---

### Script 2 — Universe Reduction via Correlation Filtering

Reduces the universe to remove redundant, highly collinear assets. Uses a graph-based approach:

- Builds an adjacency matrix where assets are connected if |ρ_{ij}| > 0.5.
- Interprets the matrix as a graph (via `igraph`) and runs a greedy independent set algorithm.
- The resulting filtered set preserves sector and style diversity while eliminating correlation clustering.

Computes eigenvalue diagnostics (condition number κ = λ_max / λ_min) to assess numerical stability of the sample covariance, motivating the regularized estimators in Scripts 3–5.

---

### Script 3 — Single Index Model (SIM) Covariance

Defines an equal-weighted market factor from filtered returns:

```
R_{m,t} = (1/N) * Σ R_{i,t}
```

Regresses each asset on the market factor to extract α_i, β_i, and idiosyncratic variance σ²_{ε_i}. Constructs the SIM covariance matrix:

```
Cov(R_i, R_j) = β_i * β_j * σ²_m        (i ≠ j)
Var(R_i)      = β²_i * σ²_m + σ²_{ε_i}  (i = j)

Σ_SIM = β β^T σ²_m + diag(σ²_ε)
```

Heatmaps show how systematic (factor) risk generates broad covariance blocks while idiosyncratic risk remains diagonal.

---

### Script 4 — Constant Correlation Model

Assumes all pairwise correlations equal the average off-diagonal correlation (≈ 0.2125). Constructs the regularized covariance matrix:

```
Σ_const = D * C_const * D
```

where D is the diagonal matrix of sample standard deviations and C_const has 1s on the diagonal and the average correlation everywhere else. This reduces sampling noise while preserving individual volatility estimates.

---

### Script 5 — Shrinkage Covariance

Implements linear shrinkage toward the diagonal (variances only):

```
Σ^shrink = (1 - γ) * Σ + γ * D,   γ = 0.4
```

Off-diagonal covariances are scaled down by a factor of 0.6, tempering extreme co-movement estimates and improving out-of-sample stability.

---

### Script 6 — Markowitz Efficient Frontier & Tangency Portfolios

Uses Lagrange multiplier–based mean–variance optimization to construct the efficient frontier across a grid of target returns. For each covariance estimator (sample, SIM, constant-correlation, shrinkage), computes a tangency portfolio:

```
w^tang ∝ Σ⁻¹ (μ - r_f * 1)
```

Weights are normalized to sum to 1. Outputs include the efficient frontier, individual asset risk–return coordinates, and four tangency portfolios for visual comparison.

---

### Script 7 — Black–Litterman Framework

**Benchmark construction:** Averages the four tangency portfolio weights to form a consensus benchmark:

```
w^bench_i = (1/4) * Σ_{methods} w_{i, method}
```

**Equilibrium returns:** Derives the implied equilibrium excess returns that make the benchmark mean–variance optimal, using a risk-aversion parameter δ calibrated from the benchmark's own Sharpe ratio:

```
δ = (E[R_bench] - r_f) / Var(R_bench)
π_excess = δ * Σ * w_bench
```

**Investor views (P, Q, Ω):** Encodes a view on the highest-return asset, asserting its future excess return aligns with its historical sample excess return. View uncertainty is quantified in Ω.

**Posterior returns and covariance:**

```
μ^BL_excess = [(τΣ)⁻¹ + P^T Ω⁻¹ P]⁻¹ [(τΣ)⁻¹ π_excess + P^T Ω⁻¹ Q]

Σ^BL = Σ + [(τΣ)⁻¹ + P^T Ω⁻¹ P]⁻¹
```

**BL optimal portfolio:** Computes a tangency portfolio using posterior returns and covariance, then constructs a BL efficient frontier and visualizes benchmark vs. optimal positioning.

---

## Key Outputs

| Output | Description |
|---|---|
| `efficient_frontier_with_optimal_portfolios.png` | Markowitz frontier with 4 tangency portfolios overlaid |
| `bl_efficient_frontier_with_optimal.png` | BL frontier, benchmark, and BL optimal portfolio |
| `optimal_portfolios_weights.csv` | Asset weights for all four tangency portfolios |
| `bl_weights_benchmark_and_optimal.csv` | Side-by-side benchmark vs. BL optimal weights |
| `bl_portfolio_summary.csv` | Return, volatility, and Sharpe ratio for BL optimal portfolio |

---

## How to Run

Execute the scripts sequentially in R. Each script saves intermediate `.rds` objects to `Outputs/Objects/` which are consumed by later scripts.

```r
source("Scripts/Stats Project Script 1.R")  # Data prep & full Σ
source("Scripts/Stats Project Script 2.R")  # Correlation filter
source("Scripts/Stats Project Script 3.R")  # SIM covariance
source("Scripts/Stats Project Script 4.R")  # Constant-correlation Σ
source("Scripts/Stats Project Script 5.R")  # Shrinkage Σ
source("Scripts/Stats Project Script 6.R")  # Efficient frontier
source("Scripts/Stats Project Script 7.R")  # Black–Litterman
```

**Required R packages:** `xts`, `zoo`, `ggplot2`, `igraph`, `readxl`, `openxlsx`, `dplyr`, `tidyr`, `reshape2`

---

## Conceptual Takeaways

1. **Structured universe reduction** — Correlation filtering trims 138 names to a numerically stable, lower-redundancy subset while preserving sector breadth.
2. **Covariance model comparison** — Sample, SIM, constant-correlation, and shrinkage estimators each make different bias–variance trade-offs; the portfolio implications are meaningfully different.
3. **Efficient frontier as a research canvas** — Tangency portfolios under different Σ reveal how risk model choice shifts optimal allocations and how diversification dominates individual securities.
4. **Black–Litterman as a Bayesian bridge** — Rather than trusting raw historical means, BL anchors to a market-consistent benchmark prior and blends in views with explicit uncertainty, producing more intuitive and stable portfolios.
5. **Extensible framework** — The pipeline can be updated with new price data and revised views without changing the underlying model structure.
