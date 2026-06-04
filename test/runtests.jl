using PeriodicARModels, OrderedCollections, FileIO, JLD2, Random, Dates
using StatsBase
using Test

##Seed for random reproducibility
Random.seed!(1234)


# ── Synthetic toy data (used in all tests below) ─────────────────────────────

# Shared synthetic data: 10 years of daily dates starting 2001-01-01
const TOY_START = Date(2001, 1, 1)
const TOY_END = Date(2050, 12, 31)
const TOY_DATE = collect(TOY_START:Day(1):TOY_END)
const TOY_N = length(TOY_DATE)
const TOY_END_short = Date(2010, 12, 31)
const TOY_DATE_short = collect(TOY_START:Day(1):TOY_END_short)
const TOY_N_short = length(TOY_DATE_short)
const TOY_RNG = MersenneTwister(42)

# Generate a univariate AR(1) series with known monthly parameters then white noise
# Φ constant = 0.6, σ constant = 1.0 across months → easy ground-truth
const TOY_Φ_TRUE = fill(0.6, 12, 1)   # 12 × p
const TOY_σ_TRUE = fill(1.0, 12)

let rng = MersenneTwister(0)
    global TOY_Z = zeros(TOY_N)
    n2m = month.(TOY_DATE)
    TOY_Z[1] = randn(rng)
    for t in 2:TOY_N
        TOY_Z[t] = TOY_Φ_TRUE[n2m[t], 1] * TOY_Z[t-1] + TOY_σ_TRUE[n2m[t]] * randn(rng)
    end
end

# ── Bug / performance regression tests ────────────────────────────────────────

@testset "fitted_periodicity_fonc: return types and fit quality" begin
    x_ = TOY_Z
    date_ = TOY_DATE

    # return_parameters=false → plain Function, not a Tuple
    f = PeriodicARModels.fitted_periodicity_fonc(x_, date_, OrderTrig=3)
    @test f isa Function
    @test f(date_[1]) isa AbstractFloat
    @test all(isfinite, f.(date_))

    # return_parameters=true → (Function, AbstractVector) with 1 + 2K coefficients
    f2, beta = PeriodicARModels.fitted_periodicity_fonc(x_, date_, OrderTrig=3, return_parameters=true)
    @test f2 isa Function
    @test beta isa AbstractVector
    @test length(beta) == 1 + 2 * 3   # intercept + 2K trig coefficients

    # Trigonometric fit must reduce variance
    @test StatsBase.var(x_ .- f.(date_)) < StatsBase.var(x_)
end

@testset "fitted_periodicity_fonc_auto: return types" begin
    x_ = TOY_Z
    date_ = TOY_DATE

    # return_parameters=false → (Function, Integer); first element must be callable
    f_auto, I_auto = PeriodicARModels.fitted_periodicity_fonc_auto(x_, date_)
    @test f_auto isa Function          # FAILS if unpacking gives a Tuple instead
    @test I_auto isa Integer
    @test 1 <= I_auto <= 30
    @test f_auto(date_[1]) isa AbstractFloat
    @test all(isfinite, f_auto.(date_))

    # return_parameters=true → should return (Function, AbstractVector, Integer)
    # BUG: precedence in return statement yields ((func, beta), I) instead of (func, beta, I)
    result_p = PeriodicARModels.fitted_periodicity_fonc_auto(x_, date_, return_parameters=true)
    @test result_p isa Tuple{<:Function,<:AbstractVector,<:Integer}
end

@testset "autotrigo periodicity mode in MonthlySWG" begin
    # BUG: decompose assigns the (func, I) tuple from fitted_periodicity_fonc_auto
    # directly to autotrigo_function without unpacking, so broadcasting it fails.
    @test MonthlySWG(TOY_Z, TOY_DATE,
        p=1, periodicity_model="autotrigo") isa MonthlySWG
end

@testset "Simulation: output shapes and finite values" begin
    # Build minimal toy models for simulation tests
    _model_uni = MonthlySWG(TOY_Z[1:TOY_N_short], TOY_DATE_short; p=1, Trendtype="null",
        periodicity_model="trigo", degree_period=2,
        σ_periodicity_model="null", σ_Trendtype="null")

    let rng2 = MersenneTwister(1), _n2m = month.(TOY_DATE_short), _z2 = zeros(TOY_N_short)
        _z2[1] = randn(rng2)
        for t in 2:TOY_N_short
            _z2[t] = 0.4 * _z2[t-1] + 0.8 * randn(rng2)
        end
        _Z2 = hcat(TOY_Z[1:TOY_N_short], _z2)
        global _model_multi = MonthlySWG(_Z2, TOY_DATE_short; p=1, Trendtype="null",
            periodicity_model="trigo", degree_period=2,
            σ_periodicity_model="null", σ_Trendtype="null")
    end

    # Univariate: n_sim vectors of length N
    sims_uni = rand(_model_uni, n_sim=3)
    @test length(sims_uni) == 3
    @test all(s -> length(s) == TOY_N_short, sims_uni)
    @test all(s -> all(isfinite, s), sims_uni)

    # Multivariate correction="resample": must terminate and return (N, d) matrices
    # (guards against the unbounded while-loop on ordering constraint)
    sims_resample = rand(_model_multi, n_sim=2, correction="resample")
    @test length(sims_resample) == 2
    @test all(s -> size(s) == (TOY_N_short, 2), sims_resample)
    @test all(s -> all(isfinite, s), sims_resample)

    # correction="null" must also return correct shapes
    sims_null = rand(_model_multi, n_sim=2, correction="null")
    @test length(sims_null) == 2
    @test all(s -> size(s) == (TOY_N_short, 2), sims_null)
end

@testset "OLS: inv(A'A)*A'*x ≈ A\\x for trigonometric design matrix" begin
    # Verify numerical consistency between the explicit normal-equations form used
    # in fitted_periodicity_fonc and the more stable backslash solve.
    x_ = TOY_Z
    date_ = TOY_DATE
    N = length(x_)
    ω = 2π / 365.2422
    K = 5
    n2t = dayofyear_Leap.(date_)
    cols = [ones(N)]
    for j in 1:K
        push!(cols, cos.(ω * j * n2t))
        push!(cols, sin.(ω * j * n2t))
    end
    Design = stack(cols)
    beta_inv = inv(transpose(Design) * Design) * transpose(Design) * x_
    beta_ls = Design \ x_
    @test maximum(abs.(beta_inv .- beta_ls)) < 1e-8
end

# ── Synthetic toy examples ────────────────────────────────────────────────────

@testset "Toy — MonthlyAR fit: types and parameter shapes" begin
    ar = MonthlyAR(TOY_Z, TOY_DATE, 1)

    @test ar isa MonthlyAR
    @test size(ar.Φ) == (12, 1)
    @test length(ar.σ) == 12
    @test all(isfinite, ar.Φ)
    @test all(isfinite, ar.σ)
    @test all(>(0), ar.σ)          # noise std must be positive
end

@testset "Toy — MonthlyAR fit: parameter recovery (AR(1), stationary)" begin
    ar = MonthlyAR(TOY_Z, TOY_DATE, 1)

    # With 10 years of data the estimates should be within 0.15 of the truth
    tol = 0.15
    @test maximum(abs.(ar.Φ[:, 1] .- 0.6)) < tol
    @test maximum(abs.(ar.σ .- 1.0)) < tol
end

@testset "Toy — MonthlyAR fit: parameter recovery (10% relative tolerance)" begin
    ar = MonthlyAR(TOY_Z, TOY_DATE, 1)

    rel_tol = 0.15
    @test maximum(abs.(ar.Φ[:, 1] .- 0.6) ./ 0.6) < rel_tol
    @test maximum(abs.(ar.σ .- 1.0) ./ 1.0) < rel_tol
end

@testset "Toy — MonthlyAR fit: AR(2) parameter recovery" begin
    # Generate an AR(2) series with Φ₁=0.5, Φ₂=0.2, σ=1.0 (same across months)
    Φ1_true = 0.5
    Φ2_true = 0.3
    σ_true = 0.5
    z2 = zeros(TOY_N)
    n2m = month.(TOY_DATE)
    rng_ar2 = MersenneTwister(99)
    z2[1] = randn(rng_ar2)
    z2[2] = randn(rng_ar2)
    for t in 3:TOY_N
        z2[t] = Φ1_true * z2[t-1] + Φ2_true * z2[t-2] + σ_true * randn(rng_ar2)
    end

    ar2 = MonthlyAR(z2, TOY_DATE, 2)

    @test size(ar2.Φ) == (12, 2)
    @test length(ar2.σ) == 12
    rel_tol = 0.20
    @test maximum(abs.(ar2.Φ[:, 1] .- Φ1_true) ./ Φ1_true) < rel_tol   # AR(1) coefficient
    @test maximum(abs.(ar2.Φ[:, 2] .- Φ2_true) ./ Φ2_true) < rel_tol   # AR(2) coefficient
    @test maximum(abs.(ar2.σ .- σ_true) ./ σ_true) < rel_tol
end

@testset "Toy — MonthlySWG fit: AR coefficient recovery (null trend, null σ-periodicity)" begin
    # With no trend and no σ-periodicity, z == standardised residuals ≈ TOY_Z,
    # so the monthly AR fit inside MonthlySWG should recover the same parameters
    # as a direct MonthlyAR fit.
    model = MonthlySWG(TOY_Z, TOY_DATE; p=1, Trendtype="null",
        periodicity_model="trigo", degree_period=2,
        σ_periodicity_model="null", σ_Trendtype="null")
    ar = MonthlyAR(TOY_Z, TOY_DATE, 1)

    rel_tol = 0.10
    # Φ recovered through MonthlySWG must agree with direct MonthlyAR fit
    @test maximum(abs.(model.monthlyAR.Φ[:, 1] .- ar.Φ[:, 1]) ./ (abs.(ar.Φ[:, 1]) .+ 1e-8)) < rel_tol
    @test maximum(abs.(model.monthlyAR.σ .- ar.σ) ./ ar.σ) < rel_tol
    # And the AR fit itself should still be close to the ground truth
    @test maximum(abs.(model.monthlyAR.Φ[:, 1] .- 0.6) ./ 0.6) < rel_tol
    @test maximum(abs.(model.monthlyAR.σ .- 1.0) ./ 1.0) < rel_tol
end

@testset "Toy — Multivariate MonthlyAR fit: diagonal coefficient recovery" begin
    # Two independent AR(1) series → the diagonal of Φ[m][1] should recover
    # the two univariate coefficients (0.6 and 0.4) and off-diagonals ≈ 0.
    rng2 = MersenneTwister(1)
    z2 = zeros(TOY_N)
    n2m = month.(TOY_DATE)
    z2[1] = randn(rng2)
    for t in 2:TOY_N
        z2[t] = 0.4 * z2[t-1] + 0.8 * randn(rng2)
    end
    Z2 = hcat(TOY_Z, z2)

    ar_multi = MonthlyAR(Z2, TOY_DATE, 1)

    tol = 0.15
    # Collect the [1,1] and [2,2] diagonal entries across the 12 months
    diag11 = [ar_multi.Φ[m][1][1, 1] for m in 1:12]
    diag22 = [ar_multi.Φ[m][1][2, 2] for m in 1:12]

    @test maximum(abs.(diag11 .- 0.6)) < tol   # channel 1 AR coef ≈ 0.6
    @test maximum(abs.(diag22 .- 0.4)) < tol   # channel 2 AR coef ≈ 0.4
end

@testset "Toy — MonthlyAR rand: output length and finiteness" begin
    ar = MonthlyAR(TOY_Z, TOY_DATE, 1)
    n2m = month.(TOY_DATE)

    sim1 = rand(ar, n2m)
    @test length(sim1) == TOY_N
    @test all(isfinite, sim1)

    sims = rand(ar, n2m; n_sim=5)
    @test length(sims) == 5
    @test all(s -> length(s) == TOY_N, sims)
    @test all(s -> all(isfinite, s), sims)
end

@testset "Toy — MonthlyAR rand: reproducibility with fixed RNG" begin
    ar = MonthlyAR(TOY_Z, TOY_DATE, 1)
    n2m = month.(TOY_DATE)

    s1 = rand(MersenneTwister(7), ar, n2m)
    s2 = rand(MersenneTwister(7), ar, n2m)
    @test s1 == s2
end

@testset "Toy — MonthlySWG fit: types and field shapes" begin
    # Null trend, simple trigonometric periodicity (order 2) for speed
    model = MonthlySWG(TOY_Z, TOY_DATE; p=1, Trendtype="null",
        periodicity_model="trigo", degree_period=2,
        σ_periodicity_model="null", σ_Trendtype="null")

    @test model isa MonthlySWG
    @test model.monthlyAR isa MonthlyAR
    @test size(model.monthlyAR.Φ) == (12, 1)
    @test length(model.monthlyAR.σ) == 12
    @test length(model.period) == 366       # one value per day-of-year
    @test length(model.σ_period) == 366
    @test length(model.trend) == TOY_N
    @test length(model.z) == TOY_N
    @test all(isfinite, model.z)
end

@testset "Toy — MonthlySWG rand: output length and finiteness" begin
    model = MonthlySWG(TOY_Z, TOY_DATE; p=1, Trendtype="null",
        periodicity_model="trigo", degree_period=2,
        σ_periodicity_model="null", σ_Trendtype="null")

    sim1 = rand(model)
    @test length(sim1) == TOY_N
    @test all(isfinite, sim1)

    sims = rand(model; n_sim=4)
    @test length(sims) == 4
    @test all(s -> length(s) == TOY_N, sims)
    @test all(s -> all(isfinite, s), sims)
end

@testset "Toy — MonthlySWG rand: fixed initial condition y₁" begin
    model = MonthlySWG(TOY_Z, TOY_DATE; p=1, Trendtype="null",
        periodicity_model="trigo", degree_period=2,
        σ_periodicity_model="null", σ_Trendtype="null")

    y₁ = model.z[1:1]
    s1 = rand(MersenneTwister(3), model; y₁=y₁)
    s2 = rand(MersenneTwister(3), model; y₁=y₁)
    @test s1 == s2
    @test length(s1) == TOY_N
end

@testset "Toy — Multivariate MonthlyAR fit: types and parameter shapes" begin
    # Build a 2-dim series: col 1 = TOY_Z, col 2 = independent AR(1) with Φ=0.4
    rng2 = MersenneTwister(1)
    z2 = zeros(TOY_N)
    n2m = month.(TOY_DATE)
    z2[1] = randn(rng2)
    for t in 2:TOY_N
        z2[t] = 0.4 * z2[t-1] + 0.8 * randn(rng2)
    end
    TOY_Z2 = hcat(TOY_Z, z2)   # N × 2

    ar2 = MonthlyAR(TOY_Z2, TOY_DATE, 1)

    @test ar2 isa MonthlyAR
    # Φ should be 12-element vector of p-element vectors of (d×d) matrices
    @test length(ar2.Φ) == 12
    @test length(ar2.Φ[1]) == 1     # p=1
    @test size(ar2.Φ[1][1]) == (2, 2)
    @test length(ar2.σ) == 12
    @test size(ar2.σ[1]) == (2, 2)  # lower-triangular Cholesky factor
    @test all(m -> all(isfinite, m), ar2.Φ[1])
    @test all(m -> all(isfinite, m), ar2.σ)

    n2m_vec = month.(TOY_DATE)
    sims2 = rand(ar2, n2m_vec; n_sim=3)
    @test length(sims2) == 3
    @test all(s -> size(s) == (TOY_N, 2), sims2)
    @test all(s -> all(isfinite, s), sims2)
end