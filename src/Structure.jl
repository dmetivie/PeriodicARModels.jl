using FileIO, JLD2, ConcreteStructs

abstract type AR_SWG end

@concrete struct MonthlyAR <: AR_SWG
    Φ                    # AbstractArray
    σ                    # AbstractVector
end


@concrete struct MonthlySWG <: AR_SWG
    monthlyAR
    trend                # AbstractArray
    period               # AbstractArray
    σ_trend              # AbstractArray
    σ_period             # AbstractArray
    date_vec             # AbstractArray
    z                    # AbstractArray
end

# @concrete struct Multi_MonthlyAR <: AR_SWG
#     Φ                    # AbstractArray
#     σ                    # AbstractArray
# end


# @concrete struct Multi_MonthlySWG <: AR_SWG
#     multi_MonthlyAR
#     trend                # AbstractArray
#     period               # AbstractArray
#     σ_trend              # AbstractArray
#     σ_period             # AbstractArray
#     date_vec             # AbstractArray
#     y₁                   # AbstractArray
#     z                    # AbstractArray
# end

include("utils.jl")
include("Periodicity.jl")
include("table_reader.jl")
include("Estimation.jl")
include("Simulation.jl")
include("Trend.jl")
include("Multi_AR_Estimation.jl")
# include("Plotting.jl")



# function fit_simpleAR(x, date_vec, p, periodicity_model::String, degree_period::Integer)
#     if periodicity_model == "trigo"
#         trigo_function = fitted_periodicity_fonc(x, date_vec, OrderTrig=degree_period)
#         periodicity = trigo_function.(date_vec)
#         nspart = trigo_function.(Date(0):Date(1)-Day(1))
#     elseif periodicity_model == "smooth"
#         smooth_function = fitted_smooth_periodicity_fonc(x, date_vec, OrderDiff=degree_period)
#         periodicity = smooth_function.(date_vec)
#         nspart = smooth_function.(Date(0):Date(1)-Day(1))
#     elseif periodicity_model == "autotrigo"
#         autotrigo_function = fitted_periodicity_fonc_stepwise(x, date_vec, MaxOrder=degree_period)
#         periodicity = autotrigo_function.(date_vec)
#         nspart = autotrigo_function.(Date(0):Date(1)-Day(1))
#     end
#     y = x - periodicity
#     Φ, σ = LL_AR_Estimation(y, p)
#     return SimpleAR(Φ, σ, nspart, y[1:p])
# end

# function fit_simpleAR(x, date_vec, p, periodicity_model::String="trigo")
#     if periodicity_model == "trigo"
#         return fit_simpleAR(x, date_vec, p, periodicity_model, 5)
#     elseif periodicity_model == "smooth"
#         return fit_simpleAR(x, date_vec, p, periodicity_model, 9)
#     elseif periodicity_model == "mean"
#         nspart = mean.(GatherYearScenario(x, date_vec))
#         y = x - nspart[dayofyear_Leap.(date_vec)]
#         Φ, σ = LL_AR_Estimation(y, p)
#         return SimpleAR(Φ, σ, nspart, y[1:p])
#     end
# end





# fit_simpleAR(x, date_vec, p=2, degree_period::Integer=5) = fit_simpleAR(x, date_vec, p, "trigo", degree_period)

# series = extract_series("TX_STAID000031.txt", plot=false)
# x, date_vec = (series[!, 2], series.DATE)
# myAR = fit_simpleAR(x, date_vec, 1)

# inverse_dayofyear_Leap(n) = Date(0) + Day(n - 1)

ismatrix(M) = false
ismatrix(M::AbstractMatrix) = true
ismulti(model::MonthlySWG) = ismatrix(model.period)

function random_init_cond(Φ, σ, t, rng::Random.AbstractRNG)
    if ismatrix(σ[1])
        p, d = length(Φ[1]), size(σ[1])[2]
        y₁ = stack([σ[t] * randn(rng, d) for _ in 1:p], dims=1)
    else
        y₁ = [rand(rng, Normal(0, σ[t])) for _ in 1:(size(Φ)[2])]
    end
    return y₁
end

function Base.rand(rng::Random.AbstractRNG, model::MonthlySWG, n2t::AbstractVector; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1)
    if ismatrix(model.period)
        period = model.period[dayofyear_Leap.(model.date_vec), :]
        σ_period = model.σ_period[dayofyear_Leap.(model.date_vec), :]
    else
        period = model.period[dayofyear_Leap.(model.date_vec)]
        σ_period = model.σ_period[dayofyear_Leap.(model.date_vec)]
    end
    nspart = model.trend .+ period
    σ_nspart = model.σ_trend .* σ_period
    if isnothing(y₁)
        if n_sim == 1
            return SimulateScenario(random_init_cond(model.monthlyAR.Φ, model.monthlyAR.σ, n2t[1], rng), n2t, model.monthlyAR.Φ, model.monthlyAR.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res)
        else
            if return_res
                return unzip([SimulateScenario(random_init_cond(model.monthlyAR.Φ, model.monthlyAR.σ, n2t[1], rng), n2t, model.monthlyAR.Φ, model.monthlyAR.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res) for _ in 1:n_sim])
            else
                return [SimulateScenario(random_init_cond(model.monthlyAR.Φ, model.monthlyAR.σ, n2t[1], rng), n2t, model.monthlyAR.Φ, model.monthlyAR.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res) for _ in 1:n_sim]
            end
        end
    else
        if n_sim == 1
            return SimulateScenario(y₁, n2t, model.monthlyAR.Φ, model.monthlyAR.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res)
        else
            if return_res
                return unzip([SimulateScenario(y₁, n2t, model.monthlyAR.Φ, model.monthlyAR.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res) for _ in 1:n_sim])
            else
                return [SimulateScenario(y₁, n2t, model.monthlyAR.Φ, model.monthlyAR.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res) for _ in 1:n_sim]
            end
        end
    end
end
function Base.rand(rng::Random.AbstractRNG, model::MonthlyAR, n2t::AbstractVector; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1)
    if isnothing(y₁)
        if n_sim == 1
            return (nspart == 0 && σ_nspart == 1) ? SimulateScenario(random_init_cond(model.Φ, model.σ, n2t[1], rng), n2t, model.Φ, model.σ, rng=rng, correction=correction) : SimulateScenario(random_init_cond(model.Φ, model.σ, n2t[1], rng), n2t, model.Φ, model.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res)
        else
            return (nspart == 0 && σ_nspart == 1) ? [SimulateScenario(random_init_cond(model.Φ, model.σ, n2t[1], rng), n2t, model.Φ, model.σ, rng=rng, correction=correction) for _ in 1:n_sim] : [SimulateScenario(random_init_cond(model.Φ, model.σ, n2t[1], rng), n2t, model.Φ, model.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res) for _ in 1:n_sim]
        end
    else
        if n_sim == 1
            return (nspart == 0 && σ_nspart == 1) ? SimulateScenario(y₁, n2t, model.Φ, model.σ, rng=rng, correction=correction) : SimulateScenario(y₁, n2t, model.Φ, model.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res)
        else
            return (nspart == 0 && σ_nspart == 1) ? [SimulateScenario(y₁, n2t, model.Φ, model.σ, rng=rng, correction=correction) for _ in 1:n_sim] : [SimulateScenario(y₁, n2t, model.Φ, model.σ, nspart, σ_nspart, rng=rng, correction=correction, return_res=return_res) for _ in 1:n_sim]
        end
    end
end
function Base.rand(rng::Random.AbstractRNG, model::AR_SWG, date_vec::AbstractVector{Date}; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1)
    return rand(rng, model, month.(date_vec), n_sim=n_sim, y₁=y₁, return_res=return_res, correction=correction, nspart=nspart, σ_nspart=σ_nspart)
end
Base.rand(model::AR_SWG, n2t::AbstractVector; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1) = rand(Random.default_rng(), model, n2t, n_sim=n_sim, y₁=y₁, correction=correction, return_res=return_res, nspart=nspart, σ_nspart=σ_nspart)
Base.rand(model::AR_SWG, date_vec::AbstractVector{Date}; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1) = rand(Random.default_rng(), model, month.(date_vec), n_sim=n_sim, y₁=y₁, correction=correction, return_res=return_res, nspart=nspart, σ_nspart=σ_nspart)
Base.rand(rng::Random.AbstractRNG, model::MonthlySWG; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1) = rand(rng, model, month.(model.date_vec), n_sim=n_sim, y₁=y₁, correction=correction, return_res=return_res, nspart=nspart, σ_nspart=σ_nspart)
Base.rand(model::AR_SWG; n_sim::Integer=1, y₁=nothing, correction="resample", return_res=false, nspart=0, σ_nspart=1) = rand(Random.default_rng(), model, month.(model.date_vec), n_sim=n_sim, y₁=y₁, correction=correction, return_res=return_res, nspart=nspart, σ_nspart=σ_nspart)

# rand(myAR, date_vec[1]:date_vec[end], 100, y₁=0.1)


function fit_ARMonthlyParameters(y, date_vec, p, method_, Nb_try=0)
    method_ == "monthlyLL" ? nothing : Monthly_temp = MonthlySeparateX(y, date_vec)
    if method_ == "mean"
        Monthly_Estimators = MonthlyEstimation(Monthly_temp, p) #Monthly_Estimators[i][j][k][l] i-> month, j-> year, k-> 1 for [Φ_1,Φ_2,...], 2 for σ, l -> index of the parameter (Φⱼ) of year if k=1 
        Monthly_Estimators2 = [[[year_[1]; year_[2]] for year_ in Month] |> stack for Month in Monthly_Estimators]
        meanparam = mulmean.(eachrow.(Monthly_Estimators2)) |> stack
        return meanparam[1:2, :]', meanparam[3, :]
    elseif method_ == "median"
        Monthly_Estimators = MonthlyEstimation(Monthly_temp, p) #Monthly_Estimators[i][j][k][l] i-> month, j-> year, k-> 1 for [Φ_1,Φ_2,...], 2 for σ, l -> index of the parameter (Φⱼ) of year if k=1 
        Monthly_Estimators2 = [[[year_[1]; year_[2]] for year_ in Month] |> stack for Month in Monthly_Estimators]
        medianparam = mulmedian.(eachrow.(Monthly_Estimators2)) |> stack
        return medianparam[1:2, :]', medianparam[3, :]
    elseif method_ == "concat"
        return MonthlyConcatanatedEstimation(Monthly_temp, p) #Φ[i][j] : i-> month, j-> index of the parameter (Φⱼ) or (σ)
    elseif method_ == "sumLL"
        return MonthlyEstimationSumLL(Monthly_temp, p)
    elseif method_ == "monthlyLL"
        return LL_AR_Estimation_monthly(y, date_vec, p, Nb_try=Nb_try)
    end
end


function MonthlyAR(z::AbstractVector, date_vec::AbstractVector{Date},
    p::Integer=1;
    method_="monthlyLL",
    Nb_try=0)

    isnothing(method_) ? method_ = "monthlyLL" : nothing
    Φ, σ = fit_ARMonthlyParameters(z, date_vec, p, method_, Nb_try)

    return MonthlyAR(Φ, σ)
end
function MonthlyAR(z::AbstractMatrix, date_vec::AbstractVector{Date},
    p::Integer=1;
    method_="monthly",
    Nb_try=0)

    if method_ == "monthly" || isnothing(method_)
        Φ, Σ = ParseMonthlyParameter(LL_Multi_AR_Estimation_monthly(z, date_vec, p), size(z)[2])
    else
        nothing #Take Φ Σ daily (e.g 366-length list of matrix for Σ)
    end

    return MonthlyAR(Φ, Σ)
end


function MonthlySWG(x, date_vec;
    p::Integer=1,
    method_=nothing,
    periodicity_model::String="trigo",
    degree_period::Integer=0,
    Trendtype="LOESS",
    trendparam=nothing,
    σ_periodicity_model::String="trigo",
    σ_degree_period::Integer=0,
    σ_Trendtype="LOESS",
    σ_trendparam=nothing,
    Nb_try=0)

    z, trend, period, σ_trend, σ_period = decompose(x, date_vec, periodicity_model, degree_period, Trendtype, trendparam, σ_periodicity_model, σ_degree_period, σ_Trendtype, σ_trendparam)

    AR_model = MonthlyAR(z, date_vec, p, method_=method_, Nb_try=Nb_try)

    return MonthlySWG(AR_model, trend, period, σ_trend, σ_period, date_vec, z)
end


defaultparam = Dict([("LOESS", 0.08), ("polynomial", 1), ("null", 1)])
defaultorder = Dict([("trigo", 5), ("smooth", 9), ("autotrigo", 50), ("stepwise_trigo", 50)])
function decompose(x, date_vec,
    periodicity_model::String="trigo",
    degree_period::Integer=0,
    Trendtype="LOESS",
    trendparam=nothing,
    σ_periodicity_model::String="trigo",
    σ_degree_period::Integer=0,
    σ_Trendtype="LOESS",
    σ_trendparam=nothing)

    isnothing(trendparam) ? trendparam = defaultparam[Trendtype] : nothing
    degree_period == 0 ? degree_period = defaultorder[periodicity_model] : nothing
    isnothing(σ_trendparam) ? σ_trendparam = defaultparam[σ_Trendtype] : nothing
    σ_periodicity_model != "null" ? σ_degree_period == 0 ? σ_degree_period = defaultorder[σ_periodicity_model] : nothing : nothing

    z, trend_mat, period_mat, σ_trend_mat, σ_period_mat = AbstractVector[], AbstractVector[], AbstractVector[], AbstractVector[], AbstractVector[]

    for x_ in eachcol(x)
        if Trendtype == "LOESS"
            trend = LOESS(x_, trendparam)
        elseif Trendtype == "polynomial"
            trend = PolyTrendFunc(x_, trendparam).(eachindex(x))
        else
            trend = zero(x_)
        end
        y_ = x_ - trend
        push!(trend_mat, trend)

        if periodicity_model == "trigo"
            trigo_function = fitted_periodicity_fonc(y_, date_vec, OrderTrig=degree_period)
            periodicity, period = trigo_function.(date_vec), trigo_function.(Date(0):(Date(1)-Day(1)))
        elseif periodicity_model == "smooth"
            smooth_function = fitted_smooth_periodicity_fonc(y_, date_vec, OrderDiff=degree_period)
            periodicity, period = smooth_function.(date_vec), smooth_function.(Date(0):(Date(1)-Day(1)))
        elseif periodicity_model == "autotrigo"
            autotrigo_function = fitted_periodicity_fonc_auto(y_, date_vec, MaxOrder=degree_period)
            periodicity, period = autotrigo_function.(date_vec), autotrigo_function.(Date(0):(Date(1)-Day(1)))
        end
        z_ = y_ - periodicity
        push!(period_mat, period)


        if σ_Trendtype != "null"

            if σ_Trendtype == "LOESS"
                σ_trend_sq = LOESS(z_ .^ 2, σ_trendparam)
            elseif σ_Trendtype == "polynomial"
                σ_trend_sq = PolyTrendFunc(z_ .^ 2, σ_trendparam).(eachindex(x))
            end
            σ_trend = σ_trend_sq .^ 0.5
        else
            σ_trend = ones(length(z_))
        end
        z_ = z_ ./ σ_trend
        push!(σ_trend_mat, σ_trend)


        if σ_periodicity_model == "trigo"
            trigo_function = fitted_periodicity_fonc(z_ .^ 2, date_vec, OrderTrig=σ_degree_period)
            σ_periodicity, σ_period = trigo_function.(date_vec) .^ 0.5, trigo_function.(Date(0):(Date(1)-Day(1))) .^ 0.5
        elseif σ_periodicity_model == "smooth"
            smooth_function = fitted_smooth_periodicity_fonc(z_ .^ 2, date_vec, OrderDiff=σ_degree_period)
            σ_periodicity, σ_period = smooth_function.(date_vec) .^ 0.5, smooth_function.(Date(0):(Date(1)-Day(1))) .^ 0.5
        elseif σ_periodicity_model == "autotrigo"
            autotrigo_function = fitted_periodicity_fonc_auto(z_ .^ 2, date_vec, MaxOrder=σ_degree_period)
            σ_periodicity, σ_period = autotrigo_function.(date_vec) .^ 0.5, autotrigo_function.(Date(0):(Date(1)-Day(1))) .^ 0.5
        else
            σ_periodicity, σ_period = ones(length(z_)), ones(366)
        end
        z_ = z_ ./ σ_periodicity
        push!(σ_period_mat, σ_period)
        push!(z, z_)
    end
    return length(eachcol(x)) == 1 ? (z[1], trend_mat[1], period_mat[1], σ_trend_mat[1], σ_period_mat[1]) : stack.((z, trend_mat, period_mat, σ_trend_mat, σ_period_mat))
end

save_model(model, title="model.jld2") = save(title, "model", model)
load_model(file) = load(file)["model"]

@concrete struct CaracteristicsSeries
    avg_day
    max_day
    df_month
    p
    k
end

function CaracteristicsSeries(df, p=0, k=0)
    Days_list = GatherYearScenario(df[!,2], df.DATE)
    avg_day = mean.(Days_list)
    max_day = maximum.(Days_list)
    df_month = @chain df begin
        @transform(:TEMP = df[!,2]) #Give a common name for TX, TN, etc...
        @transform(:MONTH = month.(:DATE)) #add month column
        @by(:MONTH, :MONTHLY_MEAN = mean(:TEMP), :MONTHLY_STD = std(:TEMP), :MONTHLY_MAX = maximum(:TEMP)) # grouby MONTH + takes the mean/std in each category 
    end
    return CaracteristicsSeries(avg_day, max_day, df_month, p, k)
end
CaracteristicsSeries(df, typesdata::AbstractVector{String}, p=0, k=0) = CaracteristicsSeries.([select(df, "DATE", typedata) for typedata in typesdata], p, k)