# include("utils.jl")
# include("../presentation/presutils.jl")
# include("ACF_PACF.jl")

using CairoMakie, OrderedCollections, StatsBase

##### PLOTTING #####
"""
    MiddleMonth(year::Int)

Return the indexes of the middle of each month of the input year. For now it's useful only to display month labels in plots.
"""
MiddleMonth(year::Integer) = (cumsum(DaysPerMonth(year)) .+ cumsum([0; DaysPerMonth(year)])[1:12]) ./ 2

"""
    PlotYearCurves(curvesvec::AbstractVector, labelvec::AbstractVector, title::String="", bands::AbstractVector=[], colorbands::AbstractVector=[])

Plot the annual series in curvesvec. Their length must be 365 or 366. You can also plots bands : each band must be a tuple or a vector of two series, the first one corresponding to the bottom of the band.
You can choose the color of each band in the vector colorbands.
"""
function PlotYearCurvesAxes!(fig, curvesvec::AbstractVector, title::String="", bands::AbstractVector=Tuple[], colorbands::AbstractVector=Tuple[]; colors::AbstractVector=String[], ylabel=true)
    length(curvesvec) != 0 ? length(curvesvec[1]) == 1 ? curvesvec = [curvesvec] : nothing : nothing  #We test if curvesvec is one series or a vector of series
    n_days = length(curvesvec) != 0 ? length(curvesvec[1]) : length(bands[1][1])
    ReferenceYear = n_days == 366 ? 0 : 1
    ax2 = Axis(fig, xticks=(MiddleMonth(ReferenceYear), Month_vec_low),
        ygridvisible=false,
        yticksvisible=false,
        yticklabelsvisible=false,
        xgridvisible=false,
        xticksvisible=false,
        # xticklabelspace=5.0,
        xticklabelsize=18)
    # xticklabelrotation=45.0)
    ax2.limits = ([0; n_days], nothing)
    ax = Axis(fig, xticks=[0; cumsum(DaysPerMonth(ReferenceYear))], xticklabelsvisible=false)
    pltbands = Plot[]
    for (band_, colorband) in zip(bands, colorbands)
        push!(pltbands, CairoMakie.band!(ax, 1:n_days, band_[1], band_[2]; color=colorband))
    end
    pltlines = Plot[]
    if length(colors) != 0
        for (vec, color_) in zip(curvesvec, colors)
            push!(pltlines, CairoMakie.lines!(ax, 1:n_days, vec, color=color_))
        end
    else
        for vec in curvesvec
            push!(pltlines, CairoMakie.lines!(ax, 1:n_days, vec))
        end
    end
    pltvec = [pltlines; pltbands]
    ax.title = title
    ax.titlesize = 19
    # ax.xlabel = "Day"
    # ax.xticklabelspace = 20.0
    if ylabel
        ax.ylabel = "Temperature (°C)"
        ax.ylabelsize = 18
    end

    ax.xticks = [0; cumsum(DaysPerMonth(ReferenceYear))]
    ax.limits = ([0; n_days], nothing)
    ax.yticklabelsize = 18
    return pltvec
end

function PlotYearCurves(curvesvec::AbstractVector, labelvec::AbstractVector, title::String="", bands::AbstractVector=Tuple[], colorbands::AbstractVector=Tuple[]; colors::AbstractVector=String[])
    fig = Figure(size=(900, 750))
    pltvec = PlotYearCurvesAxes!(fig[1:2, 1:2], curvesvec, title, bands, colorbands; colors=colors)
    Legend(fig[3, 1:2], pltvec, labelvec)
    return fig
end









"""
    PlotParameters(Parameters_vec::AbstractVector)

With a list of parameters estimated with the three methods studied in this project, this function return a figure object with all the parameters estimated.
You can add the real parameters values, but it is optional.
The input must be like this : [[Φ1_month_vec,Φ1_month_concat,Φ1_month_sumLL,Φ1_month_MLL,Φ1_true_param],[Φ2_month_vec,Φ2_month_concat,Φ2_month_sumLL,Φ2_month_MLL,Φ2_true_param],...,[σ_month_vec,σ_month_concat,σ_month_sumLL,σ_month_MLL,σ_true_param]]
"""
function PlotParametersMLL(Parameters_vec::AbstractVector, lines_::Bool=false)
    fig = Figure(size=(800, 600 * length(Parameters_vec)))
    month_vec, month_concat, month_sumLL, month_MLL, true_param = 0, 0, 0, 0, 0
    for (j, Parameters) in enumerate(Parameters_vec)
        try
            month_vec, month_concat, month_sumLL, month_MLL, true_param = Parameters
        catch BoundsError
            month_vec, month_concat, month_sumLL, month_MLL = Parameters
            true_param = nothing
        end
        if j == length(Parameters_vec) #i.e Parameter = σ
            for i in 1:12
                month_vec[i] = month_vec[i][month_vec[i].>1e-2] #I remove the values estimated close to 0.
            end
        end
        month_vec[1] = month_vec[1][abs.(month_vec[1]).<1e5] #To remove problematical values (exploded values)
        ax, plt1 = CairoMakie.boxplot(fig[1+3(j-1):2+3(j-1), 1:2], fill(1, length(month_vec[1])), month_vec[1]; width=0.3, color="orange")
        for i in 2:12
            month_vec[i] = month_vec[i][abs.(month_vec[i]).<1e5]
            CairoMakie.boxplot!(ax, fill(i, length(month_vec[i])), month_vec[i]; width=0.3, color="orange")
        end
        lines_ ? pltl = lines!(ax, collect(1:12), j == length(Parameters_vec) ? mean.(month_vec) : median.(month_vec); color="blue") : nothing
        plt2 = isnothing(true_param) ? nothing : CairoMakie.scatter!(ax, collect(1:12), true_param; color="Blue", markersize=15)
        I_concat = findall(abs.(month_concat) .< 1e4) #To remove problematical values (exploded values)
        plt3 = CairoMakie.scatter!(ax, I_concat .+ 0.15, month_concat[I_concat]; color="red", marker=:utriangle, markersize=12.5)
        I_sumLL = findall(abs.(month_sumLL) .< 1e4) #To remove problematical values (exploded values)
        plt4 = CairoMakie.scatter!(ax, I_sumLL .- 0.15, month_sumLL[I_sumLL]; color="green", marker=:dtriangle, markersize=12.5)
        plt5 = CairoMakie.scatter!(ax, collect(1:12), month_MLL; color="Purple", markersize=12.5)
        str = j == length(Parameters_vec) ? "σ" : "Φ$(j)"
        ax.title = isnothing(true_param) ? "Estimated $(str) with 3 methods" : "Real $(str) vs estimated $(str) with 3 methods"
        ax.xticks = (1:12, Month_vec)
        ax.xticklabelrotation = 45.0
        j == length(Parameters_vec) ? ax.ylabel = "Temperature (°C)" : nothing
        strline = j == length(Parameters_vec) ? "Mean of $(str) estimated on each year and month" : "Median of $(str) estimated on each year and month"
        if isnothing(true_param)
            if lines_
                Legend(fig[3+3(j-1), 1:2], [plt1, pltl, plt3, plt4, plt5], ["Boxplots of $(str) estimated on each year and month", strline, "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods", "Estimated $(str) with monthly likelihood"])
            else
                Legend(fig[3+3(j-1), 1:2], [plt1, plt3, plt4, plt5], ["Boxplots of $(str) estimated on each year and month", "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods", "Estimated $(str) with monthly likelihood"])
            end
        else
            if lines_
                Legend(fig[3+3(j-1), 1:2], [plt2, plt1, pltl, plt3, plt4, plt5], ["Real $(str)", "Boxplots of $(str) estimated on each year and month", strline, "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods", "Estimated $(str) with monthly likelihood"])
            else
                Legend(fig[3+3(j-1), 1:2], [plt2, plt1, plt3, plt4, plt5], ["Real $(str)", "Boxplots of $(str) estimated on each year and month", "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods", "Estimated $(str) with monthly likelihood"])
            end
        end
    end
    return fig
end
function PlotParameters(Parameters_vec::AbstractVector, lines_::Bool=false, MLL::Bool=false)
    if MLL
        return PlotParametersMLL(Parameters_vec, lines_)
    else
        fig = Figure(size=(800, 600 * length(Parameters_vec)))
        month_vec, month_concat, month_sumLL, true_param = 0, 0, 0, 0
        for (j, Parameters) in enumerate(Parameters_vec)
            try
                month_vec, month_concat, month_sumLL, true_param = Parameters
            catch BoundsError
                month_vec, month_concat, month_sumLL = Parameters
                true_param = nothing
            end
            if j == length(Parameters_vec) #i.e Parameter = σ
                for i in 1:12
                    month_vec[i] = month_vec[i][month_vec[i].>1e-2] #I remove the values estimated close to 0.
                end
            end
            month_vec[1] = month_vec[1][abs.(month_vec[1]).<1e5] #To remove problematical values (exploded values)
            ax, plt1 = CairoMakie.boxplot(fig[1+3(j-1):2+3(j-1), 1:2], fill(1, length(month_vec[1])), month_vec[1]; width=0.3, color="orange")
            for i in 2:12
                month_vec[i] = month_vec[i][abs.(month_vec[i]).<1e5]
                CairoMakie.boxplot!(ax, fill(i, length(month_vec[i])), month_vec[i]; width=0.3, color="orange")
            end
            lines_ ? pltl = lines!(ax, collect(1:12), j == length(Parameters_vec) ? mean.(month_vec) : median.(month_vec); color="blue") : nothing
            plt2 = isnothing(true_param) ? nothing : CairoMakie.scatter!(ax, collect(1:12), true_param; color="Blue", markersize=15)
            I_concat = findall(abs.(month_concat) .< 1e4) #To remove problematical values (exploded values)
            plt3 = CairoMakie.scatter!(ax, I_concat .+ 0.15, month_concat[I_concat]; color="red", marker=:utriangle, markersize=12.5)
            I_sumLL = findall(abs.(month_sumLL) .< 1e4) #To remove problematical values (exploded values)
            plt4 = CairoMakie.scatter!(ax, I_sumLL .- 0.15, month_sumLL[I_sumLL]; color="green", marker=:dtriangle, markersize=12.5)
            str = j == length(Parameters_vec) ? "σ" : "Φ$(j)"
            ax.title = isnothing(true_param) ? "Estimated $(str) with 3 methods" : "Real $(str) vs estimated $(str) with 3 methods"
            ax.xticks = (1:12, Month_vec)
            ax.xticklabelrotation = 45.0
            j == length(Parameters_vec) ? ax.ylabel = "Temperature (°C)" : nothing
            strline = j == length(Parameters_vec) ? "Mean of $(str) estimated on each year and month" : "Median of $(str) estimated on each year and month"
            if isnothing(true_param)
                if lines_
                    Legend(fig[3+3(j-1), 1:2], [plt1, pltl, plt3, plt4], ["Boxplots of $(str) estimated on each year and month", strline, "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods"])
                else
                    Legend(fig[3+3(j-1), 1:2], [plt1, plt3, plt4], ["Boxplots of $(str) estimated on each year and month", "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods"])
                end
            else
                if lines_
                    Legend(fig[3+3(j-1), 1:2], [plt2, plt1, pltl, plt3, plt4], ["Real $(str)", "Boxplots of $(str) estimated on each year and month", strline, "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods"])
                else
                    Legend(fig[3+3(j-1), 1:2], [plt2, plt1, plt3, plt4], ["Real $(str)", "Boxplots of $(str) estimated on each year and month", "Estimated $(str) with months concatanated", "Estimated $(str) with sum of likelihoods"])
                end
            end
        end
        return fig
    end
end



"""
    PlotMonthlyStats(RealStats::AbstractVector,SimulatedStats::AbstractMatrix,Stats::String)

Plot the monthly statistics in RealStats, the range of the monthly stats from simulations in SimulatedStats (row=month,column=simulations) and the quantile interval (0.25,0.75) of the monthly stats from these simulations.
"""
function PlotMonthlyStatsAx!(fig, RealStats::AbstractVector, SimulatedStats::AbstractMatrix, Stats::String, comment=nothing; ylabel=true)
    ax = Axis(fig)
    ax.title = isnothing(comment) ? "Real monthly $(Stats) vs range of simulated monthly $(Stats)" : "$(comment) Real monthly $(Stats) vs\n range of simulated monthly $(Stats)"
    ax.xticks = (1:12, Month_vec_low)
    ax.xticklabelsize = 18
    ax.ylabelsize = 18
    ax.yticklabelsize = 18
    ax.titlesize = 19
    # ax.xticklabelrotation = 45.0
    ylabel ? ax.ylabel = "Temperature (°C)" : nothing
    pltvec = Plot[]
    bands = [(minimum.(eachrow(SimulatedStats)), maximum.(eachrow(SimulatedStats))),
        (quantile.(eachrow(SimulatedStats), 0.25), quantile.(eachrow(SimulatedStats), 0.75))
    ]
    colorbands = [("#009bff", 0.2), ("#009bff", 0.5)]
    for (band_, colorband) in zip(bands, colorbands)
        push!(pltvec, CairoMakie.band!(ax, 1:12, band_[1], band_[2]; color=colorband))
    end
    push!(pltvec, CairoMakie.scatter!(ax, RealStats, color="Orange"))
    return pltvec
end

function PlotMonthlyStats(RealStats::AbstractVector, SimulatedStats::AbstractMatrix, Stats::String, comment=nothing)
    fig = Figure(size=(900, 750))
    pltvec = PlotMonthlyStatsAx!(fig[1:2, 1:2], RealStats, SimulatedStats, Stats, comment)
    Legend(fig[3, 1:2], pltvec, ["Range of simulated monthly $(Stats)", "Simulated monthly $(Stats) quantile interval, p ∈ [0.25,0.75]", "Real monthly $(Stats)"])
    return fig
end


"""
    WrapPlotMonthlyStats(df_month::DataFrame,sample_::AbstractVector,sample_timeline::AbstractVector{Date})

A wrapper for PlotMonthlyStats, where the real monthly statistics are in df_month (row=month), sample_ is a vector containing the simulations and sample_timeline a vector of the timeline of the simulations.
Return three plots, for the mean, the standard deviation and the maximum.
"""
function WrapPlotMonthlyStats(df_month::DataFrame, sample_::AbstractVector, sample_timeline::AbstractVector{Date}, comment=nothing)
    idx_m = [findall(month.(sample_timeline) .== m) for m in 1:12]
    mean_ts = [[mean(ts[idx_m[m]]) for m in 1:12] for ts in sample_] |> stack
    std_ts = [[std(ts[idx_m[m]]) for m in 1:12] for ts in sample_] |> stack
    max_ts = [[maximum(ts[idx_m[m]]) for m in 1:12] for ts in sample_] |> stack
    return PlotMonthlyStats(df_month.MONTHLY_MEAN, mean_ts, "mean", comment),
    PlotMonthlyStats(df_month.MONTHLY_STD, std_ts, "standard deviation", comment),
    PlotMonthlyStats(df_month.MONTHLY_MAX, max_ts, "maximum", comment)
end

# "Monthly $(Stats) parameters"
function PlotMonthlyStatsAx(subfig, Stats_vec::AbstractVector, Stats::String; unit="", title="Estimation mensuelle de $(Stats)", color="purple")
    ax = Axis(subfig)
    ax.title = title
    ax.titlesize = 16
    ax.xticks = (1:12, Month_vec2)
    ax.ylabel = unit != "" ? Stats * " (" * unit * ")" : Stats
    ax.xticklabelrotation = 45.0
    scatterlines!(ax, 1:12, Stats_vec, color=color)
    return ax
end

function PlotMonthlyStatsAx(subfig, Stats_vec_mat::Matrix, Stats::String; unit="", title="Monthly $(Stats) parameters", colors=nothing)
    ax = Axis(subfig)
    ax.title = title
    ax.titlesize = 22
    ax.xticks = (1:12, Month_vec2)
    ax.ylabel = unit != "" ? Stats * " (" * unit * ")" : Stats
    if isnothing(colors)
        for Stats_vec in eachcol(Stats_vec_mat)
            plt = scatterlines!(ax, 1:12, Stats_vec)
        end
    else
        for (Stats_vec, color) in zip(eachcol(Stats_vec_mat), colors)
            plt = scatterlines!(ax, 1:12, Stats_vec, color=color)
        end
    end
    return ax
end




"""
    Plot the monthly parameters in the vector MonthlyParams (MonthlyParams = [[Φ1Jan,Φ1Feb...],...[σJan,σFeb]])
"""
function PlotMonthlyparams(MonthlyParams::AbstractVector{Matrix}, legend)
    fig = Figure(size=(100 + 420 * length(MonthlyParams), 400))
    for (i, MonthlyParam) in enumerate(MonthlyParams)
        if i == length(MonthlyParams)
            _, plt = PlotMonthlyStatsAx(fig[1, i], MonthlyParam, "σ", unit="C°")
        else
            PlotMonthlyStatsAx(fig[1, i], MonthlyParam, "Φ$(i)")
        end
    end
    Legend(figvec[2, eachindex(MonthlyParams)], plt, legend)
    return fig
end

function PlotMonthlyparams(MonthlyParams)
    fig = Figure(size=(100 + 420 * length(MonthlyParams), 400))
    for (i, MonthlyParam) in enumerate(MonthlyParams)
        if i == length(MonthlyParams)
            PlotMonthlyStatsAx(fig[1, i], MonthlyParam, "σ", unit="C°")
        else
            PlotMonthlyStatsAx(fig[1, i], MonthlyParam, "Φ$(i)")
        end
    end
    return fig
end

istickable10year(date_) = month(date_) == 1 && day(date_) == 1 && (year(date_) % 10 == 0)
strfunc(date_) = "$(year(date_))"

function PlotTrend_simple(Model)
    n_days = length(Model.date_vec)
    # xlimits = [-Int(floor(0.02*n_days)), n_days + Int(floor(0.05*n_days))]
    ticksindexes = findall(istickable10year, Model.date_vec)
    xticklabel = strfunc.(Model.date_vec[ticksindexes])

    fig = Figure(size=(1200, 400))

    ax1, plt = lines(fig[1, 1], Model.trend, color="green")
    # lines!(ax1, zeros(length(Model.trend)), color="black")
    ax1.limits = ([1, n_days], nothing)
    ax1.xticks = (ticksindexes, xticklabel)
    ax1.ylabel = "Temperature (°C)"
    ax1.ylabelsize = 25
    ax1.yticklabelsize = 20
    ax1.xticklabelsize = 20
    ax1.title = "Tendance additive Mₜ"
    ax1.titlesize = 26

    ax2, plt = lines(fig[1, 2], Model.σ_trend, color="green")
    ax2.limits = ([1, n_days], nothing)
    ax2.xticks = (ticksindexes, xticklabel)
    ax2.yticklabelsize = 20
    ax2.xticklabelsize = 20
    ax2.title = "Tendance multiplicative σMₜ"
    ax2.titlesize = 26
    return fig

end
function PlotTrend_multi(Model, typedata2="TX")
    n_days = length(Model.date_vec)
    # xlimits = [-Int(floor(0.02*n_days)), n_days + Int(floor(0.05*n_days))]
    ticksindexes = findall(istickable10year, Model.date_vec)
    xticklabel = strfunc.(Model.date_vec[ticksindexes])

    fig = Figure(size=(1200, 1000))

    ax1, plt = lines(fig[1, 1], Model.trend[:, 1] .- minimum(Model.trend[:, 1]), color="green")
    lines!(ax1, zeros(length(Model.trend)), color="black")
    ax1.limits = ([1, n_days], nothing)
    ax1.xticks = (ticksindexes, xticklabel)
    ax1.ylabel = "Temperature (°C)"
    ax1.ylabelsize = 25
    ax1.yticklabelsize = 20
    ax1.xticklabelsize = 20
    ax1.title = "Additive Trend of TN"
    ax1.titlesize = 26

    ax2, plt = lines(fig[1, 2], Model.σ_trend[:, 1], color="green")
    ax2.limits = ([1, n_days], nothing)
    ax2.xticks = (ticksindexes, xticklabel)
    ax2.yticklabelsize = 20
    ax2.xticklabelsize = 20
    ax2.title = "Multiplicative Trend of TN"
    ax2.titlesize = 26

    ax3, plt = lines(fig[2, 1], Model.trend[:, 2] .- minimum(Model.trend[:, 2]), color="green")
    lines!(ax3, zeros(length(Model.trend)), color="black")
    ax3.limits = ([1, n_days], nothing)
    ax3.xticks = (ticksindexes, xticklabel)
    ax3.ylabel = "Temperature (°C)"
    ax3.ylabelsize = 25
    ax3.yticklabelsize = 20
    ax3.xticklabelsize = 20
    ax3.title = "Additive Trend of $(typedata2)"
    ax3.titlesize = 26

    ax4, plt = lines(fig[2, 2], Model.σ_trend[:, 2], color="green")
    ax4.limits = ([1, n_days], nothing)
    ax4.xticks = (ticksindexes, xticklabel)
    ax4.yticklabelsize = 20
    ax4.xticklabelsize = 20
    ax4.title = "Multiplicative Trend of $(typedata2)"
    ax4.titlesize = 26

    return fig
end
PlotTrend(Model::MonthlySWG, typedata2="TX") = ismulti(Model) ? PlotTrend_multi(Model, typedata2) : PlotTrend_simple(Model)


function PlotσTrend(Model::MonthlySWG, typedata2="TG")
    if !ismulti(Model)
        error("The input model must be a Multi SWG model.")
    end
    n_days = length(Model.date_vec)
    # xlimits = [-Int(floor(0.02*n_days)), n_days + Int(floor(0.05*n_days))]
    ticksindexes = findall(istickable10year, Model.date_vec)
    xticklabel = strfunc.(Model.date_vec[ticksindexes])

    fig = Figure(size=(1200, 500))


    ax2, plt = lines(fig[1, 1], Model.σ_trend[:, 1], color="green")
    ax2.limits = ([1, n_days], nothing)
    ax2.xticks = (ticksindexes, xticklabel)
    ax2.yticklabelsize = 20
    ax2.xticklabelsize = 20
    ax2.title = "Multiplicative Trend of TN"
    ax2.titlesize = 26


    ax4, plt = lines(fig[1, 2], Model.σ_trend[:, 2], color="green")
    ax4.limits = ([1, n_days], nothing)
    ax4.xticks = (ticksindexes, xticklabel)
    ax4.yticklabelsize = 20
    ax4.xticklabelsize = 20
    ax4.title = "Multiplicative Trend of $(typedata2)"
    ax4.titlesize = 26

    return fig
end




function PlotSeasonnality_simple(Model::MonthlySWG)
    fig = Figure(size=(1200, 400))

    PlotYearCurvesAxes!(fig[1, 1], [Model.period], "Saisonnalité additive Sₜ", colors=["orange"])
    PlotYearCurvesAxes!(fig[1, 2], [Model.σ_period], "Saisonnalité multiplicative σSₜ", colors=["orange"], ylabel=false)

    return fig

end

function PlotSeasonnality_multi(Model::MonthlySWG, typedata2="TX")
    fig = Figure(size=(1200, 1000))

    PlotYearCurvesAxes!(fig[1, 1], [Model.period[:, 1]], "Additive seasonality of TN", colors=["orange"])
    PlotYearCurvesAxes!(fig[1, 2], [Model.σ_period[:, 1]], "Multiplicative seasonality of TN", colors=["orange"])
    PlotYearCurvesAxes!(fig[2, 1], [Model.period[:, 2]], "Additive seasonality of $(typedata2)", colors=["orange"])
    PlotYearCurvesAxes!(fig[2, 2], [Model.σ_period[:, 2]], "Multiplicative seasonality of $(typedata2)", colors=["orange"])

    return fig

end
PlotSeasonnality(Model::MonthlySWG, typedata2="TX") = ismulti(Model) ? PlotSeasonnality_multi(Model, typedata2) : PlotSeasonnality_simple(Model)


function Sample_diagnostic(sample_, date_vec, period, avg_day, max_day, df_month; format_="vertical", size=format_ == "vertical" ? (1200, 1900) : (1600, 900))
    year_sample = GatherYearScenarios(sample_, date_vec)

    idx_m = [findall(month.(date_vec) .== m) for m in 1:12]
    mean_ts = [[mean(ts[idx_m[m]]) for m in 1:12] for ts in sample_] |> stack
    std_ts = [[std(ts[idx_m[m]]) for m in 1:12] for ts in sample_] |> stack
    max_ts = [[maximum(ts[idx_m[m]]) for m in 1:12] for ts in sample_] |> stack

    fig = Figure(size=size)
    legendsize = 14

    if format_ == "horizontal"
        figvec = [(fig[1:2, 1:3], fig[3, 2]),
            (fig[1:2, 4:6], fig[3, 5]),
            (fig[1:2, 7:9], fig[3, 8]),
            (fig[4:5, 1:3], fig[6, 2]),
            (fig[4:5, 4:6], fig[6, 5]),
            (fig[4:5, 7:9], fig[6, 8])]
    else
        figvec = [(fig[1:2, 1:3], fig[3, 2]),
            (fig[1:2, 4:6], fig[3, 5]),
            (fig[4:5, 1:3], fig[6, 2]),
            (fig[4:5, 4:6], fig[6, 5]),
            (fig[7:8, 1:3], fig[9, 2]),
            (fig[7:8, 4:6], fig[9, 5])]
    end

    ylabelBoolVec = format_ == "horizontal" ? [true, false, false, true, false, false] : [true, false, true, false, true, false]

    plt1 = PlotYearCurvesAxes!(figvec[1][1], [period, mean.(year_sample)], "A. Average daily temperature during a year", ylabel=ylabelBoolVec[1])
    Legend(figvec[1][2], plt1, ["Periodicity estimation", "Mean simulated temperatures"], labelsize=legendsize)

    plt2 = PlotYearCurvesAxes!(figvec[2][1], [period, avg_day, max_day],
        "B. Simulated vs recorded daily temperatures during a year",
        [(minimum.(year_sample), maximum.(year_sample)), (quantile.(year_sample, 0.25), quantile.(year_sample, 0.75))],
        [("#009bff", 0.2), ("#009bff", 0.5)],
        colors=["blue", "orange", "red"],
        ylabel=ylabelBoolVec[2])
    Legend(figvec[2][2], plt2, ["Periodicity estimation", "Average recorded temperatures", "Maximum recorded temperatures", "Simulated temperatures range", "Simulated temperatures quantile interval, p ∈ [0.25,0.75]"], labelsize=legendsize)

    plt3 = PlotYearCurvesAxes!(figvec[3][1], [maximum.(year_sample) .- minimum.(year_sample), quantile.(year_sample, 0.75) .- quantile.(year_sample, 0.25)],
        "C. Simulated temperatures interquartile range",
        ylabel=ylabelBoolVec[3])
    Legend(figvec[3][2], plt3, ["Simulated temperatures range", "Simulated temperatures interquartile range, p ∈ [0.25,0.75]"], labelsize=legendsize)

    plt4 = PlotMonthlyStatsAx!(figvec[4][1], df_month.MONTHLY_MEAN, mean_ts, "mean", "D.", ylabel=ylabelBoolVec[4])
    Legend(figvec[4][2], plt4, ["Range of simulated monthly mean", "Simulated monthly mean quantile interval, p ∈ [0.25,0.75]", "Real monthly mean"], labelsize=legendsize)

    plt5 = PlotMonthlyStatsAx!(figvec[5][1], df_month.MONTHLY_STD, std_ts, "standard deviation", "E.", ylabel=ylabelBoolVec[5])
    Legend(figvec[5][2], plt5, ["Range of simulated monthly standard deviation", "Simulated monthly standard deviation quantile interval, p ∈ [0.25,0.75]", "Real monthly standard deviation"], labelsize=legendsize)

    plt6 = PlotMonthlyStatsAx!(figvec[6][1], df_month.MONTHLY_MAX, max_ts, "maximum", "F.", ylabel=ylabelBoolVec[6])
    Legend(figvec[6][2], plt6, ["Range of simulated monthly maximum", "Simulated monthly maximum quantile interval, p ∈ [0.25,0.75]", "Real monthly maximum"], labelsize=legendsize)

    return fig
end

function PrintParams(Φ_vec, σ)
    println("Parameters : \n")
    for (i, Φ) in enumerate(Φ_vec)
        println("Φ$(i) = $(Φ)")
    end
    println("σ = $(σ) C° \n")
end


function PrintParams(Φ_vec, σ, io)
    println(io, "Parameters : \n")
    for (i, Φ) in enumerate(Φ_vec)
        println(io, "Φ$(i) = $(Φ)")
    end
    println(io, "σ = $(σ) C° \n")
end


function Sample_diagnostic_simple(sample_::Tuple, Caracteristics_Series, Model::MonthlySWG; format_="vertical", size=format_ == "vertical" ? (1200, 1900) : (1600, 900), folder=nothing, settings=nothing)
    p, k = Caracteristics_Series.p, Caracteristics_Series.k
    sample_, z_sample = sample_
    # println(z_sample[1])
    # println(sample_[1])
    fig0 = PlotTrend(Model)

    # cteParam = length(Model.monthlyAR.σ) == 1
    # cteParam ? PrintParams(Model.monthlyAR.Φ, Model.monthlyAR.σ) : 
    fig1 = PlotMonthlyparams([eachcol(Model.monthlyAR.Φ); [Model.monthlyAR.σ]])
    fig2 = Sample_diagnostic(sample_,
        Model.date_vec,
        Model.period .+ mean.(GatherYearScenario(Model.trend, Model.date_vec)),
        Caracteristics_Series.avg_day,
        Caracteristics_Series.max_day,
        Caracteristics_Series.df_month,
        format_=format_,
        size=size
    )
    fig3 = Plot_Sample_MonthlyACF(z_sample, Model.date_vec, Model.z, "p=$(p), k=$(k)")
    fig4 = Plot_Sample_MonthlyPACF(z_sample, Model.date_vec, Model.z, "p=$(p), k=$(k)")

    if !isnothing(folder)
        mkpath(folder)
        save(folder * "/Trend_k=$(k)" * ".pdf", fig0; px_per_unit=2.0)
        # cteParam ? nothing : save(folder * "/Params" * "_$(p)_$(k)" * ".pdf", fig1; px_per_unit=2.0)
        save(folder * "/Params" * "_$(p)_$(k)" * ".pdf", fig1; px_per_unit=2.0)
        save(folder * "/Sample_diagnostic" * "_$(p)_$(k)" * ".pdf", fig2; px_per_unit=2.0)
        save(folder * "/MonthlyACF" * "_$(p)_$(k)" * ".pdf", fig3; px_per_unit=2.0)
        save(folder * "/MonthlyPACF" * "_$(p)_$(k)" * ".pdf", fig4; px_per_unit=2.0)
        open(folder * "/Params_Results" * "_$(p)_$(k)" * ".txt", "a") do io
            # cteParam ? PrintParams(Model.Φ, Model.σ, io) : nothing
            if !isnothing(settings)
                println(io, "Settings :\n")
                for key in keys(settings)
                    println(io, "$(key) : $(settings[key])")
                end
                println(io, "\n\n")
            end
            # println(io, "Results :\n")
            # println(io, "Additive periodicity order : $(k)")
            # println(io, "Multiplicative periodicity order : $(Model.σ_period_order)")

            println(io, "Number of simulations for the last diagnostic : $(length(sample_)) \n \n")

        end
    end
    # return cteParam ? (fig0, fig2, fig3, fig4) : (fig0, fig1, fig2, fig3, fig4)
    return (fig0, fig1, fig2, fig3, fig4)
end
function Sample_diagnostic_simple(sample_::AbstractVector{T}, Caracteristics_Series, Model::MonthlySWG; format_="vertical", size=format_ == "vertical" ? (1200, 1900) : (1600, 900), folder=nothing, settings=nothing) where T<:AbstractVector
    nspart_ = Model.trend .+ Model.period[dayofyear_Leap.(Model.date_vec)]
    σ_nspart_ = Model.σ_trend .* Model.σ_period[dayofyear_Leap.(Model.date_vec)]
    z_sample = map(sim -> (sim .- nspart_) ./ σ_nspart_, sample_)
    return Sample_diagnostic_simple((sample_, z_sample), Caracteristics_Series, Model; format_=format_, size=size, folder=folder, settings=settings)
end


function Sample_diagnostic_multi(sample_::Tuple, Caracteristics_Series, Model::MonthlySWG; format_="vertical", size=format_ == "vertical" ? (1200, 1900) : (1600, 900), folder=nothing, settings=nothing, TG_bool=false)
    p, k = Caracteristics_Series[1].p, Caracteristics_Series[1].k #Caracteristics_Series::Vector
    sample_, z_sample = sample_

    typedata2 = TG_bool ? "TG" : "TX"

    sampleTN, z_sampleTN = [s[:, 1] for s in sample_], [s[:, 1] for s in z_sample]
    sampleTX, z_sampleTX = [s[:, 2] for s in sample_], [s[:, 2] for s in z_sample]

    fig0 = PlotTrend(Model, typedata2)
    # fig1 = PlotMonthlyparams([invert(Model.Φ); [Model.σ]])
    fig2 = Sample_diagnostic(sampleTN,
        Model.date_vec,
        Model.period[:, 1] .+ mean.(GatherYearScenarios([Model.trend[:, 1]], Model.date_vec)),
        Caracteristics_Series[1].avg_day,
        Caracteristics_Series[1].max_day,
        Caracteristics_Series[1].df_month,
        format_=format_,
        size=size
    )
    fig3 = Sample_diagnostic(sampleTX,
        Model.date_vec,
        Model.period[:, 2] .+ mean.(GatherYearScenarios([Model.trend[:, 2]], Model.date_vec)),
        Caracteristics_Series[2].avg_day,
        Caracteristics_Series[2].max_day,
        Caracteristics_Series[2].df_month,
        format_=format_,
        size=size
    )
    fig4 = Plot_Sample_MonthlyACF(z_sampleTN, Model.date_vec, Model.z[:, 1], "TN, p=$(p), k=$(k)")
    fig5 = Plot_Sample_MonthlyPACF(z_sampleTN, Model.date_vec, Model.z[:, 1], "TN, p=$(p), k=$(k)")
    fig6 = Plot_Sample_MonthlyACF(z_sampleTX, Model.date_vec, Model.z[:, 2], "$(typedata2), p=$(p), k=$(k)")
    fig7 = Plot_Sample_MonthlyPACF(z_sampleTX, Model.date_vec, Model.z[:, 2], "$(typedata2), p=$(p), k=$(k)")
    fig8 = Plot_Sample_MonthlyCC(z_sampleTN, z_sampleTX, Model.date_vec, Model.z[:, 1], Model.z[:, 2], "p=$(p), k=$(k)", typedata2)

    if !isnothing(folder)
        mkpath(folder)
        # save(folder * "/Params" * "_$(p)_$(k)" * ".pdf", fig1; px_per_unit=2.0)
        save(folder * "/Trend_k=$(k)" * ".pdf", fig0; px_per_unit=2.0)
        save(folder * "/Sample_diagnostic_TN" * "_$(p)_$(k)" * ".pdf", fig2; px_per_unit=2.0)
        save(folder * "/Sample_diagnostic_$(typedata2)" * "_$(p)_$(k)" * ".pdf", fig3; px_per_unit=2.0)
        save(folder * "/MonthlyACF_TN" * "_$(p)_$(k)" * ".pdf", fig4; px_per_unit=2.0)
        save(folder * "/MonthlyPACF_TN" * "_$(p)_$(k)" * ".pdf", fig5; px_per_unit=2.0)
        save(folder * "/MonthlyACF_$(typedata2)" * "_$(p)_$(k)" * ".pdf", fig6; px_per_unit=2.0)
        save(folder * "/MonthlyPACF_$(typedata2)" * "_$(p)_$(k)" * ".pdf", fig7; px_per_unit=2.0)
        save(folder * "/MonthlyCC" * "_$(p)_$(k)" * ".pdf", fig8; px_per_unit=2.0)
        open(folder * "/Params_Results" * "_$(p)_$(k)" * ".txt", "a") do io
            if !isnothing(settings)
                println(io, "Settings :\n")
                for key in keys(settings)
                    println(io, "$(key) : $(settings[key])")
                end
                println(io, "\n\n")
            end
            # println(io, "Results :\n")
            # println(io, "Number of scenarios with dates where TN > $(typedata2) : $(sum(TN_Grt_TX.(sample_) .> 0))")
            # println(io, "Percentage of scenarios with dates where TN > $(typedata2) : $(trunc(100 * sum(TN_Grt_TX.(sample_) .> 0)/length(sample_),digits=2)) %")
            # println(io, "Mean percentage of dates where TN > $(typedata2) : $(trunc(100*mean(TN_Grt_TX.(sample_))/length(Model.date_vec),digits=2)) %")
            # println(io, "Median percentage of dates where TN > $(typedata2) : $(trunc(100*median(TN_Grt_TX.(sample_))/length(Model.date_vec),digits=2)) %")
            # println(io, "Additive periodicity order : $(k)")
            # println(io, "Multiplicative periodicity order : $(Model.σ_period_order)")
        end
    end
    return (fig0, fig2, fig3, fig4, fig5, fig6, fig7, fig8)
end
function Sample_diagnostic_multi(sample_::AbstractVector{T}, Caracteristics_Series, Model::MonthlySWG; format_="vertical", size=format_ == "vertical" ? (1200, 1900) : (1600, 900), folder=nothing, settings=nothing, TG_bool=false) where T<:AbstractMatrix
    nspart_ = Model.trend .+ Model.period[dayofyear_Leap.(Model.date_vec), :]
    σ_nspart_ = Model.σ_trend .* Model.σ_period[dayofyear_Leap.(Model.date_vec), :]
    z_sample = map(sim -> (sim .- nspart_) ./ σ_nspart_, sample_)
    return Sample_diagnostic_multi((sample_, z_sample), Caracteristics_Series, Model; format_=format_, size=size, folder=folder, settings=settings, TG_bool=TG_bool)
end

function Sample_diagnostic(sample_, Caracteristics_Series, Model::MonthlySWG; format_="vertical", size=format_ == "vertical" ? (1200, 1900) : (1600, 900), folder=nothing, settings=nothing, TG_bool=false)
    ismulti(Model) ? Sample_diagnostic_multi(sample_, Caracteristics_Series, Model; format_=format_, size=size, folder=folder, settings=settings, TG_bool=TG_bool) :
    Sample_diagnostic_simple(sample_, Caracteristics_Series, Model; format_=format_, size=size, folder=folder, settings=settings)
end