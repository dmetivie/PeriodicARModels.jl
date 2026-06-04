using Dates, LinearAlgebra, DataInterpolations, RegularizationTools, GLM

"""
    fitted_periodicity_fonc(x::AbstractVector,return_parameters::Bool=false)

Return a trigonometric function f of period 365.25 of equation f(t) = μ + a*cos(2π*t/365.25) + b*sin((2π*t/365.25) fitted on x. 
If return_parameters=true, return a tuple with f and [μ,a,b]. Be careful : the function returned takes the same arguments as dayofyear_Leap() (Either Date of Integer and Date, see above).
"""
function fitted_periodicity_fonc(x::AbstractVector, date_vec::AbstractVector; OrderTrig::Integer=1, return_parameters::Bool=false)
    N = length(x)
    n2t = dayofyear_Leap.(date_vec)
    ω = 2π / 365.2422
    cos_nj = [cos.(ω * j * n2t) for j = 1:OrderTrig]
    sin_nj = [sin.(ω * j * n2t) for j = 1:OrderTrig]
    Design = stack([[ones(N)]; interleave2(cos_nj, sin_nj)])
    beta = Design \ x
    function func(args...)
        t = dayofyear_Leap(args...)
        IL = interleave2([cos(ω * j * t) for j = 1:OrderTrig], [sin(ω * j * t) for j = 1:OrderTrig])
        return dot(beta, [1; IL])
    end
    return return_parameters ? (func, beta) : func
end

AIC_seas(n, p, SRS) = 2p + n * (log(2π * SRS / n) + 1)
BIC_seas(n, p, SRS) = p*log(n) + n * (log(2π * SRS / n) + 1)

"""
    fitted_smooth_periodicity_fonc(x::AbstractVector, date_vec::AbstractVector, orderdiff::Integer=9)

Return a function which is the smooth regularization of the mean year of x, with the timeline date_vec. 
"""
function fitted_smooth_periodicity_fonc(x::AbstractVector, date_vec::AbstractVector; OrderDiff::Integer=9)
    f = RegularizationSmooth(mean.(GatherYearScenario(x, date_vec)), 1:366, OrderDiff)
    return date -> f(dayofyear_Leap(date))
end


function fitted_periodicity_fonc_auto(x::AbstractVector, date_vec::AbstractVector; MaxOrder::Integer=30, return_parameters::Bool=false, Verbose::Bool=false, UltraVerbose=false)
    N = length(x)
    n2t = dayofyear_Leap.(date_vec)
    ω = 2π / 365.2422
    cos_nj = [cos.(ω * j * n2t) for j = 1:MaxOrder]
    sin_nj = [sin.(ω * j * n2t) for j = 1:MaxOrder]
    Design = stack([[ones(N)]; interleave2(cos_nj, sin_nj)])
    # Compute QR once on the full design matrix and re-use columns for each order
    F = qr(Design)
    AIC_seas_vec, beta_vec = AbstractFloat[], AbstractVector[]
    for i in 1:MaxOrder
        SubDesign = Design[:, 1:(1+2i)]
        beta = SubDesign \ x
        AIC_seas_ = AIC_seas(N, 1 + 2i, sum((SubDesign * beta .- x) .^ 2))
        UltraVerbose ? println("i=$(i), AIC_seas=$(trunc(AIC_seas_,digits=2))") : nothing
        push!(AIC_seas_vec, AIC_seas_)
        push!(beta_vec, beta)
    end
    I = argmin(AIC_seas_vec)
    beta = beta_vec[I]
    function func(args...)
        t = dayofyear_Leap(args...)
        IL = interleave2([cos(ω * j * t) for j = 1:I], [sin(ω * j * t) for j = 1:I])
        return dot(beta, [1; IL])
    end
    Verbose ? println("Best AIC ($(trunc(minimum(AIC_seas_vec)))) reached for i=$(I)") : nothing
    # Bug fix: explicit parentheses to avoid return-precedence issue with the comma
    return return_parameters ? (func, beta, I) : (func, I)
end