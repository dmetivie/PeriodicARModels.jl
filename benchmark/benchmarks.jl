using BenchmarkTools
using PeriodicARModels
using Dates
using Random
using Statistics

const SUITE = BenchmarkGroup()

# ── Shared synthetic data ─────────────────────────────────────────────────────
# Same data generation as in test/runtests.jl (TOY_* constants)

const BENCH_START = Date(2001, 1, 1)
const BENCH_DATE_10Y = collect(BENCH_START:Day(1):Date(2010, 12, 31))   # ~3 650 days
const BENCH_DATE_30Y = collect(BENCH_START:Day(1):Date(2030, 12, 31))   # ~10 958 days
const N = 300

function _simulate_ar1(dates, Φ=0.6, σ=1.0; seed=0)
    rng = MersenneTwister(seed)
    n2m = month.(dates)
    z   = zeros(length(dates))
    z[1] = randn(rng)
    for t in 2:length(dates)
        z[t] = Φ * z[t-1] + σ * randn(rng)
    end
    return z
end

const Z_UNI_10Y = _simulate_ar1(BENCH_DATE_10Y)
const Z_UNI_30Y = _simulate_ar1(BENCH_DATE_30Y)

function _simulate_var1(dates, Φ=0.4, σ=0.8; seed=1)
    rng  = MersenneTwister(seed)
    z2   = zeros(length(dates))
    z2[1] = randn(rng)
    for t in 2:length(dates)
        z2[t] = Φ * z2[t-1] + σ * randn(rng)
    end
    return hcat(_simulate_ar1(dates), z2)
end

const Z_MULTI_10Y = _simulate_var1(BENCH_DATE_10Y)
const Z_MULTI_30Y = _simulate_var1(BENCH_DATE_30Y)

# ── 1. MonthlyAR fit (univariate) ─────────────────────────────────────────────
SUITE["MonthlyAR_fit"] = BenchmarkGroup()

for (label, z, dates) in [
        ("10y_p1", Z_UNI_10Y, BENCH_DATE_10Y),
        ("10y_p2", Z_UNI_10Y, BENCH_DATE_10Y),
        ("30y_p1", Z_UNI_30Y, BENCH_DATE_30Y),
    ]
    p = endswith(label, "p1") ? 1 : 2
    SUITE["MonthlyAR_fit"][label] = @benchmarkable(
        MonthlyAR($z, $dates, $p),
        evals = 1
    )
end

# ── 2. MonthlyAR rand (univariate) ────────────────────────────────────────────
SUITE["MonthlyAR_rand"] = BenchmarkGroup()

let
    ar_10y  = MonthlyAR(Z_UNI_10Y, BENCH_DATE_10Y, 1)
    n2m_10y = month.(BENCH_DATE_10Y)
    ar_30y  = MonthlyAR(Z_UNI_30Y, BENCH_DATE_30Y, 1)
    n2m_30y = month.(BENCH_DATE_30Y)

    SUITE["MonthlyAR_rand"]["single_10y"] = @benchmarkable(
        rand($ar_10y, $n2m_10y),
        evals = 1
    )
    SUITE["MonthlyAR_rand"]["10sim_10y"] = @benchmarkable(
        rand($ar_10y, $n2m_10y; n_sim = N),
        evals = 1
    )
    SUITE["MonthlyAR_rand"]["100sim_10y"] = @benchmarkable(
        rand($ar_10y, $n2m_10y; n_sim = N),
        evals = 1
    )
    SUITE["MonthlyAR_rand"]["single_30y"] = @benchmarkable(
        rand($ar_30y, $n2m_30y),
        evals = 1
    )
end

# ── 3. MonthlySWG fit (univariate) ────────────────────────────────────────────
# Vary: periodicity model and AR order p
SUITE["MonthlySWG_fit"] = BenchmarkGroup()

for (label, kw) in [
        ("null_trend_trigo2_p1", (Trendtype="null",  periodicity_model="trigo", degree_period=2, σ_Trendtype="null",  σ_periodicity_model="null", p=1)),
        ("null_trend_trigo5_p2", (Trendtype="null",  periodicity_model="trigo", degree_period=5, σ_Trendtype="null",  σ_periodicity_model="null", p=2)),
        ("loess_trigo5_p1",      (Trendtype="LOESS", periodicity_model="trigo", degree_period=5, σ_Trendtype="LOESS", σ_periodicity_model="trigo", σ_degree_period=5, p=1)),
        ("loess_trigo5_p2",      (Trendtype="LOESS", periodicity_model="trigo", degree_period=5, σ_Trendtype="LOESS", σ_periodicity_model="trigo", σ_degree_period=5, p=2)),
    ]
    SUITE["MonthlySWG_fit"][label] = @benchmarkable(
        MonthlySWG($Z_UNI_10Y, $BENCH_DATE_10Y; $kw...),
        evals = 1
    )
end

# ── 4. MonthlySWG rand (univariate) ───────────────────────────────────────────
SUITE["MonthlySWG_rand"] = BenchmarkGroup()

let
    model_simple = MonthlySWG(Z_UNI_10Y, BENCH_DATE_10Y; p=1,
        Trendtype="null", periodicity_model="trigo", degree_period=2,
        σ_Trendtype="null", σ_periodicity_model="null")
    model_full   = MonthlySWG(Z_UNI_10Y, BENCH_DATE_10Y; p=2,
        Trendtype="LOESS", periodicity_model="trigo", degree_period=5,
        σ_Trendtype="LOESS", σ_periodicity_model="trigo", σ_degree_period=5)

    SUITE["MonthlySWG_rand"]["simple_single"] = @benchmarkable(
        rand($model_simple),
        evals = 1
    )
    SUITE["MonthlySWG_rand"]["simple_100sim"] = @benchmarkable(
        rand($model_simple; n_sim = N),
        evals = 1
    )
    SUITE["MonthlySWG_rand"]["full_single"] = @benchmarkable(
        rand($model_full),
        evals = 1
    )
    SUITE["MonthlySWG_rand"]["full_100sim"] = @benchmarkable(
        rand($model_full; n_sim = N),
        evals = 1
    )
end

# ── 5. MonthlyAR fit (multivariate, d=2) ─────────────────────────────────────
SUITE["MonthlyAR_multi_fit"] = BenchmarkGroup()

for (label, z, dates) in [
        ("10y_p1", Z_MULTI_10Y, BENCH_DATE_10Y),
        ("10y_p2", Z_MULTI_10Y, BENCH_DATE_10Y),
        ("30y_p1", Z_MULTI_30Y, BENCH_DATE_30Y),
    ]
    p = endswith(label, "p1") ? 1 : 2
    SUITE["MonthlyAR_multi_fit"][label] = @benchmarkable(
        MonthlyAR($z, $dates, $p),
        evals = 1
    )
end

# ── 6. MonthlyAR rand (multivariate, d=2) ────────────────────────────────────
SUITE["MonthlyAR_multi_rand"] = BenchmarkGroup()

let
    ar2      = MonthlyAR(Z_MULTI_10Y, BENCH_DATE_10Y, 1)
    n2m_10y  = month.(BENCH_DATE_10Y)

    SUITE["MonthlyAR_multi_rand"]["single_no_correction"] = @benchmarkable(
        rand($ar2, $n2m_10y; correction = "null"),
        evals = 1
    )
    SUITE["MonthlyAR_multi_rand"]["single_resample"] = @benchmarkable(
        rand($ar2, $n2m_10y; correction = "resample"),
        evals = 1
    )
    SUITE["MonthlyAR_multi_rand"]["10sim_resample"] = @benchmarkable(
        rand($ar2, $n2m_10y; n_sim = N, correction = "resample"),
        evals = 1
    )
end

# ── 7. MonthlySWG fit (multivariate, d=2) ────────────────────────────────────
SUITE["MonthlySWG_multi_fit"] = BenchmarkGroup()

for (label, kw) in [
        ("trigo2_p1", (Trendtype="null",  periodicity_model="trigo", degree_period=2, σ_Trendtype="null", σ_periodicity_model="null", p=1)),
        ("trigo5_p1", (Trendtype="LOESS", periodicity_model="trigo", degree_period=5, σ_Trendtype="LOESS", σ_periodicity_model="trigo", σ_degree_period=5, p=1)),
    ]
    SUITE["MonthlySWG_multi_fit"][label] = @benchmarkable(
        MonthlySWG($Z_MULTI_10Y, $BENCH_DATE_10Y; $kw...),
        evals = 1
    )
end

# ── 8. MonthlySWG rand (multivariate, d=2) ────────────────────────────────────
SUITE["MonthlySWG_multi_rand"] = BenchmarkGroup()

let
    model2 = MonthlySWG(Z_MULTI_10Y, BENCH_DATE_10Y; p=1,
        Trendtype="null", periodicity_model="trigo", degree_period=2,
        σ_Trendtype="null", σ_periodicity_model="null")

    SUITE["MonthlySWG_multi_rand"]["single_no_correction"] = @benchmarkable(
        rand($model2; correction = "null"),
        evals = 1
    )
    SUITE["MonthlySWG_multi_rand"]["single_resample"] = @benchmarkable(
        rand($model2; correction = "resample"),
        evals = 1
    )
    SUITE["MonthlySWG_multi_rand"]["10sim_resample"] = @benchmarkable(
        rand($model2; n_sim = N, correction = "resample"),
        evals = 1
    )
end

a = run(SUITE)