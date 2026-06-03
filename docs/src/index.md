```@meta
CurrentModule = PeriodicARModels
```

# PeriodicARModels

Documentation for [PeriodicARModels](https://github.com/ArnaudG0649/PeriodicARModels.jl).

A package to define, fit, and simulate AR models with **periodic (monthly) parameters**, primarily designed for Stochastic Weather Generators (SWG). Both univariate and multivariate (multi-site) series are supported.

---

## Mathematical Model

### Periodic AR($p$) — Univariate

Let $x_t$ be a daily time series (e.g. temperature). The model decomposes $x_t$ as:

```math
x_t = \mu_t + \sigma_t \cdot z_t
```

where $\mu_t$ is a deterministic non-stationary part (trend + seasonal mean) and $\sigma_t$ is a time-varying scale. The residual $z_t$ follows a **Periodic AR($p$)** process:

```math
z_t = \sum_{i=1}^{p} \Phi_{m(t),\, i} \; z_{t-i} + \sigma_{m(t)} \; \varepsilon_t, \qquad \varepsilon_t \overset{\text{iid}}{\sim} \mathcal{N}(0,1)
```

where $m(t) \in \{1,\ldots,12\}$ is the **month** of day $t$. Each month has its own independent set of AR coefficients $(\Phi_{m,1}, \ldots, \Phi_{m,p})$ and noise standard deviation $\sigma_m$.

### Periodic VAR($p$) — Multivariate ($d$ sites)

For a $d$-dimensional series $\mathbf{x}_t \in \mathbb{R}^d$, the residual process becomes:

```math
\mathbf{z}_t = \sum_{i=1}^{p} \boldsymbol{\Phi}_{m(t),\, i} \; \mathbf{z}_{t-i} + \boldsymbol{\Sigma}_{m(t)} \; \boldsymbol{\varepsilon}_t, \qquad \boldsymbol{\varepsilon}_t \overset{\text{iid}}{\sim} \mathcal{N}(\mathbf{0}, I_d)
```

where $\boldsymbol{\Phi}_{m,i} \in \mathbb{R}^{d \times d}$ are matrix-valued AR coefficients and $\boldsymbol{\Sigma}_m$ is a lower-triangular Cholesky factor so that the noise covariance is $\boldsymbol{\Sigma}_m \boldsymbol{\Sigma}_m^\top$.

---

## Non-Stationary Decomposition

The non-stationary part $\mu_t = \text{trend}_t + \text{period}_t$ is estimated in two steps:

1. **Trend** — removed first, using one of:
   - `"LOESS"` (default, `span = 0.08`) — local regression smoother
   - `"polynomial"` — global polynomial fit of a given degree
   - `"null"` — no trend removal

2. **Seasonal periodicity** — fitted on the detrended series, using one of:
   - `"trigo"` (default, order 5) — trigonometric regression:
     ```math
     \text{period}(t) = \sum_{j=1}^{K} \left[ a_j \cos\!\left(\tfrac{2\pi j t}{365.2422}\right) + b_j \sin\!\left(\tfrac{2\pi j t}{365.2422}\right) \right]
     ```
   - `"smooth"` (order 9) — regularised smoothing via `RegularizationTools`
   - `"autotrigo"` — trigonometric fit with order $K$ selected automatically by AIC

The same pipeline is applied to the variance ($\sigma_t$) independently.

---

## Data Structures

### `MonthlyAR`

Stores the purely periodic AR parameters:

| Field | Type | Description |
|---|---|---|
| `Φ` | `Matrix` (12 × p) or `Vector{Vector{Matrix}}` | AR coefficients per month. Univariate: `Φ[m, i]`; multivariate: `Φ[m][i]` is a $d\times d$ matrix |
| `σ` | `Vector` (length 12) or `Vector{Matrix}` | Noise scale per month. Univariate: scalar $\sigma_m$; multivariate: lower-triangular $\boldsymbol{\Sigma}_m$ |

### `MonthlySWG`

Full Stochastic Weather Generator wrapping a `MonthlyAR` with the non-stationary components:

| Field | Description |
|---|---|
| `monthlyAR` | The fitted `MonthlyAR` on the residual series $z$ |
| `trend` | Estimated trend $\text{trend}_t$ |
| `period` | Seasonal cycle (366-length or 366×d array) |
| `σ_trend` | Trend of the variance |
| `σ_period` | Seasonal cycle of the variance |
| `date_vec` | The training date vector |
| `z` | Standardised residual series |

---

## Fitting a Model

### High-level constructors

**Univariate `MonthlySWG`** (full pipeline — decompose then fit AR):

```julia
model = MonthlySWG(x, date_vec;
    p                   = 1,           # AR order
    periodicity_model   = "trigo",     # "trigo" | "smooth" | "autotrigo"
    degree_period       = 5,           # trigonometric order K (0 = auto-default)
    Trendtype           = "LOESS",     # "LOESS" | "polynomial" | "null"
    trendparam          = nothing,     # LOESS span or polynomial degree
    σ_periodicity_model = "trigo",     # same options for variance periodicity
    σ_degree_period     = 5,
    σ_Trendtype         = "LOESS",
)
```

**`MonthlyAR`** (fit AR on an already-standardised residual series `z`):

```julia
ar = MonthlyAR(z, date_vec, p;
    method_ = "monthlyLL"   # estimation method (see below)
)
```

### Estimation methods (`method_` keyword)

| Value | Description |
|---|---|
| `"monthlyLL"` (default) | Joint MLE over all 12 months simultaneously using `Optimization.jl` + `LBFGS` with `ForwardDiff` gradients |
| `"concat"` | Concatenate all data for the same month across years, then fit one AR per month |
| `"mean"` / `"median"` | Fit AR year-by-year per month, then average/median the parameters |
| `"sumLL"` | Maximise the sum of per-year log-likelihoods per month |

### Multivariate (`d > 1`)

Pass a matrix `x` of shape `(N, d)`:

```julia
ar_multi = MonthlyAR(z_matrix, date_vec, p)  # fits monthly VAR(p) via MLE
```

---

## Simulation

Once a model is fitted, new scenarios are generated with `rand`:

```julia
# Single scenario over a date vector
scenario = rand(model, date_vec)

# n_sim scenarios
scenarios = rand(model, date_vec; n_sim = 100)

# With a fixed RNG for reproducibility
using Random
scenarios = rand(MersenneTwister(42), model, date_vec; n_sim = 50)
```

For multivariate models a `correction = "resample"` option (default) rejects draws that violate ordering constraints between dimensions (e.g. $T_N \leq T_X$).

---

## Utilities

| Function | Description |
|---|---|
| `dayofyear_Leap(d)` | Leap-year-aware day-of-year index in $[1, 366]$ |
| `GatherYearScenario(x, date_vec)` | Group a series by day-of-year into a 366-element vector of vectors |
| `RootAR(Φ)` | Roots of the AR characteristic polynomial (stationarity check: all roots outside unit circle) |
| `AIC_seas(n, p, RSS)` / `BIC_seas(n, p, RSS)` | Information criteria for selecting trigonometric periodicity order |
| `decompose(x, date_vec, ...)` | Low-level decomposition returning `(z, trend, period, σ_trend, σ_period)` |

