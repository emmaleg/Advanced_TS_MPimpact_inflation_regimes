# Advanced_TS_MPimpact_inflation_regimes
Final project of the course "Advanced Time Series: macroeconometrics" from Master In Economics 2nd year, given by Prof. Giovanni Ricco

Replicate Canova Forero 2025 working paper (last version know: April 2025)
And add some extensions if time

## Replication notes (assumptions & implementation choices)

This repository implements a Threshold-BVAR with volatility-in-mean and stochastic volatility, following Canova & Forero (paper specification and Appendix A). The goal is to replicate the posterior of the inflation threshold P* and the implied regime indicator (Figures 3 and 4), and then extend to IRFs and additional figures.

### Important: paper version used
- The code is aligned as closely as possible with the publicly available version that contains the model equations, the regime definition, the Gibbs sampling blocks, and the impact matrix structure (Appendix A.8).
- The December 2025 PDF hosted on Google Drive could not be programmatically accessed in the assistant environment; therefore some calibration details may differ from that latest version. If the Dec-2025 version uses different sample windows, priors, or lag choices, update the config files accordingly.

### Data sources and transformations
- Data are downloaded from FRED using the `fredgraph.csv?id=SERIES` endpoint (no API key required).
- Variables are constructed as:
  - Industrial production (INDPRO): YoY log growth, `log(x_t) - log(x_{t-12})`
  - PCE price index (PCEPI): YoY log inflation, `log(x_t) - log(x_{t-12})`
  - Unemployment rate (UNRATE): level, divided by 100 (decimal)
  - Fed funds rate (FEDFUNDS): level, divided by 100 (decimal)
  - Term spread: GS10 - TB3MS, both divided by 100 (decimal)
  - Money (M2SL): YoY log growth
  - Commodity price index: proxied by PPIACO (All Commodities PPI), YoY log growth
  - Stock market: SP500 (daily), aggregated to monthly by "last observation", YoY log growth
- Sample window (monthly): 1960m1–2023m6 (configurable).

### Model choices not explicitly pinned down in the visible paper text
The paper provides the structure but not all numeric calibration choices in the accessible text. For those, the repository uses standard defaults, all centralized in `config/*.m` and easy to change:

1) VAR lag length p (monthly):
- Default: p = 12 (common monthly choice).

2) Number of lags of log-volatility in the mean J:
- Default: J = 2 (small, parsimonious; changeable).

3) Threshold delay d:
- Drawn from {1,...,dmax}, with dmax = 6 (as stated in the paper).
- Prior on d is uniform.

4) Priors (numerical hyperparameters):
- Coefficients (Phi): Minnesota-like diagonal normal prior (tightness and decay set in config).
- Impact matrix parameters (alpha): Normal prior N(0, 10 I).
- Diagonal variance terms (sigma^2): independent inverse-gamma with (a0,b0) specified in config.
- Volatility process h_t = ln(lambda_t): AR(1) with priors on (mu, F, Q) using Normal/Truncated Normal/Inverse-Gamma.

5) Single-move SV update:
- The paper motivates a "single-move" update for h_t due to the nonlinear appearance of lambda_t.
- Implementation here uses a single-move Metropolis step with proposal equal to the conditional AR(1) prior at each t, and acceptance based on the local likelihood window affected by h_t in the mean (from t to t+J).
- This is faithful to the "single-move" logic but may differ from the exact bounded proposal described if additional details are present in the Dec-2025 version.

### Identification / restrictions
- The impact matrix A_i is parameterized exactly as in Appendix A.8 using a 22-parameter vector alpha_i.
- Sign/zero restrictions from Table 1 are included as an optional check but are turned OFF by default (config flag). This keeps the initial goal focused on P* and regime inference (Figures 3 and 4).

### Outputs
- Posterior draws are saved to `_results/posterior_draws.mat`.
- Figure scripts read only saved outputs and create `_results/figures/fig2.pdf`, `fig3.pdf`, `fig4.pdf`.


