using Dates, LinearAlgebra, DataInterpolations, RegularizationTools, GLM, Loess

function PolyTrendFunc(x, order, index=eachindex(x); return_parameters=false)
    Design = [index .^ i for i in 0:order] |> stack
    beta = Design \ x
    f(t) = dot(beta, [t^i for i in 0:order])
    return return_parameters ? (f, beta) : f
end

# ========== LOESS ========== #

LOESS(x, span=0.4, degree=1) = predict(loess(eachindex(x), x, span=span > 1 ? span / length(x) : span, degree=degree), eachindex(x))

#!!Savoir interpreter span