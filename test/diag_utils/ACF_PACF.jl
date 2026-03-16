using CairoMakie, Statistics, StatsBase, PeriodicARModels

include("Separation.jl")

invert_v(L::AbstractVector) = [[L[i][j] for i in eachindex(L)] for j in eachindex(L[1])]

##### ACF/PACF #####
"""
    ACF_PACF(x::AbstractVector, return_data::Bool=false)

Return the graphs of the ACF and the PACF of the series x. 
If return_data = true, the function return a tuple with the ACF, the PACF and the figure object which contains the graphs.
"""
function ACF_PACF(x::AbstractVector, return_data::Bool=false)
    fig = Figure(size=(700, 800))
    autocor_ = autocor(x, 0:15)
    ax1, plt1 = barplot(fig[1, 1], 0:15, autocor_)
    ax1.title = "ACF"
    ax1.xticks = 0:15
    ax1.xgridvisible = false
    ax1.yticks = round.(range(minimum(autocor_), maximum(autocor_), 10), digits=2)
    pacf_ = pacf(x, 1:15)
    ax2, plt2 = barplot(fig[2, 1], 1:15, pacf_)
    ax2.title = "PACF"
    ax2.xticks = 1:15
    ax2.xgridvisible = false
    ax2.yticks = round.(range(minimum(pacf_), maximum(pacf_), 10), digits=2)
    return return_data ? (autocor_, pacf_, fig) : fig
end

"""
    MonthlyACF(Monthly_temp::Vector, return_data::Bool=false)

Return the graphs of the monthly average ACF of each series inside Monthly_temp.
Monthly_temp must be a 3 level-nested vector such as an output of MonthlySeparateX.
If return_data = true, the function return a tuple with the ACF of each months inside Monthly_temp and the figure object.
"""
function MonthlyACF(Monthly_temp::AbstractVector, return_data::Bool=false)
    fig = Figure(size=(800, 600))
    supertitle = Label(fig[1, 1:4], "Monthly average ACF", fontsize=20)
    ax_vec = Axis[]
    min_y = 0
    if return_data
        autocor_data = AbstractVector[]
    end
    for i in 1:12
        autocor_vec = [autocor(Monthly_temp[i][j], 0:10) for j in eachindex(Monthly_temp[i])] #One element per year
        m_autocor = mean(autocor_vec)
        ax, plot = barplot(fig[((i-1)÷4)+2, (i-1)%4+1], 0:10, m_autocor)
        ax.title = Month_vec[i]
        min_y = min(minimum(m_autocor), min_y)
        ax_vec = [ax_vec; [ax]]
        if return_data
            append!(autocor_data, [autocor_vec])
        end
    end
    for ax in ax_vec
        ax.limits = (nothing, [min_y - 0.05, 1.05])
        ax.xgridvisible = false
        ax.xticks = 0:10
    end
    return return_data ? (autocor_data, fig) : fig
end
MonthlyACF(x::AbstractVector{T}, sample_timeline::AbstractVector, return_data=false) where T<:AbstractFloat = MonthlyACF(MonthlySeparateX(x, sample_timeline), return_data)


"""
    MonthlyPACF(Monthly_temp, return_data::Bool=false)
    
Return the graphs of the monthly average PACF of each series inside Monthly_temp.
Monthly_temp must be a 3 level-nested vector such as an output of MonthlySeparateX.
If return_data = true, the function return a tuple with the PACF of each months inside Monthly_temp and the figure object.
"""
function MonthlyPACF(Monthly_temp, return_data::Bool=false)
    fig = Figure(size=(800, 600))
    supertitle = Label(fig[1, 1:4], "Monthly average PACF", fontsize=20)
    ax_vec = Axis[]
    max_y, min_y = 0, 0
    if return_data
        pacf_data = AbstractVector[]
    end
    for i in 1:12
        pacf_vec = [pacf(Monthly_temp[i][j], 1:10) for j in eachindex(Monthly_temp[i])] #One element per year
        m_pacf = mean(pacf_vec)
        ax, plot = barplot(fig[((i-1)÷4)+2, (i-1)%4+1], 1:10, m_pacf)
        ax.title = Month_vec[i]
        max_y = max(maximum(m_pacf), max_y)
        min_y = min(minimum(m_pacf), min_y)
        ax_vec = [ax_vec; [ax]]
        if return_data
            append!(pacf_data, [pacf_vec])
        end
    end
    for ax in ax_vec
        ax.limits = (nothing, [min_y - 0.05, max_y + 0.05])
        ax.xgridvisible = false
        ax.xticks = 1:10
    end
    return return_data ? (pacf_data, fig) : fig
end
MonthlyPACF(x::AbstractVector{T}, sample_timeline::AbstractVector, return_data=false) where T<:AbstractFloat = MonthlyPACF(MonthlySeparateX(x, sample_timeline), return_data)

###Monthly ACF samples###
"""
    MatrixMonthlyACF(Monthly_temp)

On a nested list of monthly series (like an output of MonthlySeparateX()), return a matrix of the mean ACF for each month.
The 12 rows represent each month and the index of the column represents the order of the ACF (between 1 and 10). 
"""
function MatrixMonthlyACF(Monthly_temp)
    acf_vec = [mean([autocor(Monthly_temp[1][j], 1:10) for j in eachindex(Monthly_temp[1])])]
    for i in 2:12
        push!(acf_vec, mean([autocor(Monthly_temp[i][j], 1:10) for j in eachindex(Monthly_temp[i])]))
    end
    return stack(acf_vec, dims=1)
end

"""
    Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date})

Return the vector of the matrixes of the mean monthly ACF (see MatrixMonthlyACF()) of each scenario in samples, with their common timeline in sample_timeline.
"""
function Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date})
    return [MatrixMonthlyACF(MonthlySeparateX(sample_, sample_timeline)) for sample_ in samples]
end


"""
    Plot_Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp=nothing)

Plot the boxplots of the mean monthly ACF of the scenarios in samples and if asked the mean monthly ACF of the series in Monthly_temp.
"""
function Plot_Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, true_matrix=nothing, comment="")
    list_matrix = Sample_MonthlyACF(samples, sample_timeline)
    fig = Figure(size=(800, 700))
    # supertitle = Label(fig[1, 1:4], "Monthly average ACF " * comment, fontsize=20)
    ax_vec = Axis[]
    min_y, max_y = 0, 0
    for i in 1:11
        sample_acf = invert_v([matrix_[i, :] for matrix_ in list_matrix])
        ax, _ = CairoMakie.barplot(fig[((i-1)÷4)+1, (i-1)%4+1], 1:10, maximum.(sample_acf); fillto=minimum.(sample_acf), color=("#009bff", 0.5))
        scatter!(ax, 1:10, median.(sample_acf), color=("#009bff", 0.8), marker=:hline, markersize=15)
        ax.title = Month_vec[i]
        isnothing(true_matrix) ? nothing : scatter!(ax, 1:10, true_matrix[i, :], color="#e57420", marker=:hline, markersize=15)
        push!(ax_vec, ax)
        CompleteMonthlySample = [reduce(vcat, [matrix_[i, :] for matrix_ in list_matrix]); true_matrix[i, :]] #All the autocorrelations in the i-th month.
        max_y = max(maximum(CompleteMonthlySample), max_y)
        min_y = min(minimum(CompleteMonthlySample), min_y)
    end
    ##### Same part as inside the for loop but we keep plot1 and plot2 for legends.
    sample_acf = invert_v([matrix_[12, :] for matrix_ in list_matrix])
    ax, plot1 = CairoMakie.barplot(fig[((12-1)÷4)+1, (12-1)%4+1], 1:10, maximum.(sample_acf); fillto=minimum.(sample_acf), color=("#009bff", 0.5))
    plot2 = scatter!(ax, 1:10, median.(sample_acf), color=("#009bff", 0.8), marker=:hline, markersize=15)
    ax.title = Month_vec[12]
    isnothing(true_matrix) ? nothing : plot3 = scatter!(ax, 1:10, true_matrix[12, :], color="#e57420", marker=:hline, markersize=15)
    push!(ax_vec, ax)
    CompleteMonthlySample = [reduce(vcat, [matrix_[12, :] for matrix_ in list_matrix]); true_matrix[12, :]]
    max_y = max(maximum(CompleteMonthlySample), max_y)
    min_y = min(minimum(CompleteMonthlySample), min_y)
    #####
    for ax in ax_vec
        ax.limits = (nothing, [min_y - 0.15, max_y + 0.15])
        ax.xgridvisible = false
        ax.xticks = 0:10
    end
    Legend(fig[4, 1:4], [plot1, plot2, plot3], ["Range of means autocorrelations of the simulated temperatures", "Median of means autocorrelations of the simulated temperatures", "mean autocorrelation of the recorded temperatures"])
    return fig
end
Plot_Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp::AbstractVector{T}, comment="") where T<:AbstractVector = Plot_Sample_MonthlyACF(samples, sample_timeline, MatrixMonthlyACF(Monthly_temp), comment)
Plot_Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, x::AbstractVector{T}, comment="") where T<:AbstractFloat = Plot_Sample_MonthlyACF(samples, sample_timeline, MonthlySeparateX(x, sample_timeline), comment)

"""
    Plot_Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp=nothing)

Plot the boxplots of the mean monthly ACF of the scenarios in samples and if asked the mean monthly ACF of the series in Monthly_temp.
"""
function Error_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp)
    list_matrix = Sample_MonthlyACF(samples, sample_timeline)
    true_matrix = MatrixMonthlyACF(Monthly_temp)
    Error_matrix(matrix) = abs.(matrix - true_matrix)
    return mean(Error_matrix.(list_matrix))
end


###Monthly PACF samples###
"""
    MatrixMonthlyACF(Monthly_temp)

On a nested list of monthly series (like an output of MonthlySeparateX()), return a matrix of the mean PACF for each month.
The 12 rows represent each month and the index of the column represents the order of the PACF (between 1 and 10). 
"""
function MatrixMonthlyPACF(Monthly_temp)
    acf_vec = [mean([pacf(Monthly_temp[1][j], 1:10) for j in eachindex(Monthly_temp[1])])]
    for i in 2:12
        push!(acf_vec, mean([pacf(Monthly_temp[i][j], 1:10) for j in eachindex(Monthly_temp[i])]))
    end
    return stack(acf_vec, dims=1)
end
#À optimiser...
"""
    Sample_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date})

Return the vector of the matrixes of the mean monthly PACF (see MatrixMonthlyPACF()) of each scenario in samples, with their common timeline in sample_timeline.
"""
function Sample_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date})
    return [MatrixMonthlyPACF(MonthlySeparateX(sample_, sample_timeline)) for sample_ in samples]
end


"""
    Plot_Sample_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp=nothing)

Plot the boxplots of the mean monthly PACF of the scenarios in samples and if asked the mean monthly PACF of the series in Monthly_temp.
"""
function Plot_Sample_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, true_matrix=nothing, comment="")
    list_matrix = Sample_MonthlyPACF(samples, sample_timeline)
    fig = Figure(size=(800, 700))
    # supertitle = Label(fig[1, 1:4], "Monthly average PACF " * comment, fontsize=20)
    ax_vec = Axis[]
    min_y, max_y = 0, 0
    for i in 1:11
        sample_pacf = invert_v([matrix_[i, :] for matrix_ in list_matrix])
        ax, _ = CairoMakie.barplot(fig[((i-1)÷4)+1, (i-1)%4+1], 1:10, maximum.(sample_pacf); fillto=minimum.(sample_pacf), color=("#009bff", 0.5))
        scatter!(ax, 1:10, median.(sample_pacf), color=("#009bff", 0.8), marker=:hline, markersize=15)
        ax.title = Month_vec[i]
        isnothing(true_matrix) ? nothing : scatter!(ax, 1:10, true_matrix[i, :], color="#e57420", marker=:hline, markersize=15)
        push!(ax_vec, ax)
        CompleteMonthlySample = [reduce(vcat, [matrix_[i, :] for matrix_ in list_matrix]); true_matrix[i, :]] #All the pacf in the i-th month.
        max_y = max(maximum(CompleteMonthlySample), max_y)
        min_y = min(minimum(CompleteMonthlySample), min_y)
    end
    ##### Same part as inside the for loop but we keep plot1 and plot2 for legends.
    sample_pacf = invert_v([matrix_[12, :] for matrix_ in list_matrix])
    ax, plot1 = CairoMakie.barplot(fig[((12-1)÷4)+1, (12-1)%4+1], 1:10, maximum.(sample_pacf); fillto=minimum.(sample_pacf), color=("#009bff", 0.5))
    plot2 = scatter!(ax, 1:10, median.(sample_pacf), color=("#009bff", 0.8), marker=:hline, markersize=15)
    ax.title = Month_vec[12]
    isnothing(true_matrix) ? nothing : plot3 = scatter!(ax, 1:10, true_matrix[12, :], color="#e57420", marker=:hline, markersize=15)
    push!(ax_vec, ax)
    CompleteMonthlySample = [reduce(vcat, [matrix_[12, :] for matrix_ in list_matrix]); true_matrix[12, :]]
    max_y = max(maximum(CompleteMonthlySample), max_y)
    min_y = min(minimum(CompleteMonthlySample), min_y)
    #####
    for ax in ax_vec
        ax.limits = (nothing, [min_y - 0.15, max_y + 0.15])
        ax.xgridvisible = false
        ax.xticks = 0:10
    end
    Legend(fig[4, 1:4], [plot1, plot2, plot3], ["Range of means PACF of the simulated temperatures", "Median of means PACF of the simulated temperatures", "mean PACF of the recorded temperatures"])
    return fig
end
Plot_Sample_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp::AbstractVector{T}, comment="") where T<:AbstractVector = Plot_Sample_MonthlyPACF(samples, sample_timeline, MatrixMonthlyPACF(Monthly_temp), comment)
Plot_Sample_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, x::AbstractVector{T}, comment="") where T<:AbstractFloat = Plot_Sample_MonthlyPACF(samples, sample_timeline, MonthlySeparateX(x, sample_timeline), comment)


"""
    Plot_Sample_MonthlyACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp=nothing)

Plot the boxplots of the mean monthly ACF of the scenarios in samples and if asked the mean monthly ACF of the series in Monthly_temp.
"""
function Error_MonthlyPACF(samples::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_temp)
    list_matrix = Sample_MonthlyPACF(samples, sample_timeline)
    true_matrix = MatrixMonthlyPACF(Monthly_temp)
    Error_matrix(matrix) = abs.(matrix - true_matrix)
    return mean(Error_matrix.(list_matrix))
end


#####CrossCorrelation#####
CrossCorrelation(X, Y, k=4) = @views [cor(X[max(1 - i, 1):min(end - i, end)], Y[max(i + 1, 1):min(end + i, end)]) for i in -k:k]

MatrixMonthlyCC(Monthly_tempX, Monthly_tempY) = stack(
    [mean([CrossCorrelation(Monthly_tempX[i][j], Monthly_tempY[i][j]) for j in eachindex(Monthly_tempX[i])]) for i in 1:12],
    dims=1
)

function Sample_MonthlyCC(samplesTN::AbstractVector, samplesTX::AbstractVector, sample_timeline::AbstractVector{Date})
    return [MatrixMonthlyCC(MonthlySeparateX(sampleTN, sample_timeline), MonthlySeparateX(sampleTX, sample_timeline)) for (sampleTN, sampleTX) in zip(samplesTN, samplesTX)]
end


function Plot_Sample_MonthlyCC(samplesTN::AbstractVector, samplesTX::AbstractVector, sample_timeline::AbstractVector{Date}, true_matrix=nothing, comment="", typedata2="TX")
    list_matrix = Sample_MonthlyCC(samplesTN, samplesTX, sample_timeline)
    fig = Figure(size=(800, 600))
    # supertitle = Label(fig[1, 1:4], "Monthly average correlation between the residuals of TNₜ and $(typedata2)ₜ₊ₛ\n" * comment, fontsize=20)
    ax_vec = Axis[]
    min_y, max_y = 0, 0
    for i in 1:11
        sample_cc = invert_v([matrix_[i, :] for matrix_ in list_matrix])
        ax, _ = CairoMakie.barplot(fig[((i-1)÷4)+1, (i-1)%4+1], -4:4, maximum.(sample_cc); fillto=minimum.(sample_cc), color=("#009bff", 0.5))
        scatter!(ax, -4:4, median.(sample_cc), color=("#009bff", 0.8), marker=:hline, markersize=15)
        ax.title = Month_vec[i]
        i > 8 ? ax.xlabel = "s" : nothing
        isnothing(true_matrix) ? nothing : scatter!(ax, -4:4, true_matrix[i, :], color="#e57420", marker=:hline, markersize=15)
        push!(ax_vec, ax)
        CompleteMonthlySample = [reduce(vcat, [matrix_[i, :] for matrix_ in list_matrix]); true_matrix[i, :]] #All the pacf in the i-th month.
        max_y = max(maximum(CompleteMonthlySample), max_y)
        min_y = min(minimum(CompleteMonthlySample), min_y)
    end
    ##### Same part as inside the for loop but we keep plot1 and plot2 for legends.
    sample_cc = invert_v([matrix_[12, :] for matrix_ in list_matrix])
    ax, plot1 = CairoMakie.barplot(fig[((12-1)÷4)+1, (12-1)%4+1], -4:4, maximum.(sample_cc); fillto=minimum.(sample_cc), color=("#009bff", 0.5))
    plot2 = scatter!(ax, -4:4, median.(sample_cc), color=("#009bff", 0.8), marker=:hline, markersize=15)
    ax.title = Month_vec[12]
    ax.xlabel = "s"
    isnothing(true_matrix) ? nothing : plot3 = scatter!(ax, -4:4, true_matrix[12, :], color="#e57420", marker=:hline, markersize=15)
    push!(ax_vec, ax)
    CompleteMonthlySample = [reduce(vcat, [matrix_[12, :] for matrix_ in list_matrix]); true_matrix[12, :]]
    max_y = max(maximum(CompleteMonthlySample), max_y)
    min_y = min(minimum(CompleteMonthlySample), min_y)
    #####
    for ax in ax_vec
        ax.limits = (nothing, [min_y - 0.15, max_y + 0.15])
        ax.xgridvisible = false
        ax.xticks = -4:4
        ax.xticklabelsize = 11
    end
    Legend(fig[4, 1:4], [plot1, plot2, plot3], ["Range of means CC of the simulated temperatures", "Median of means CC of the simulated temperatures", "mean CC of the recorded temperatures"])
    return fig
end
Plot_Sample_MonthlyCC(samplesTN::AbstractVector, samplesTX::AbstractVector, sample_timeline::AbstractVector{Date}, Monthly_tempX::AbstractVector{T}, Monthly_tempY::AbstractVector{T}, comment="", typedata2="TX") where T<:AbstractVector = Plot_Sample_MonthlyCC(samplesTN, samplesTX, sample_timeline, MatrixMonthlyCC(Monthly_tempX, Monthly_tempY), comment, typedata2)
Plot_Sample_MonthlyCC(samplesTN::AbstractVector, samplesTX::AbstractVector, sample_timeline::AbstractVector{Date}, x::AbstractVector{T}, y::AbstractVector{T}, comment="", typedata2="TX") where T<:AbstractFloat = Plot_Sample_MonthlyCC(samplesTN, samplesTX, sample_timeline, MonthlySeparateX(x, sample_timeline), MonthlySeparateX(y, sample_timeline), comment, typedata2)
