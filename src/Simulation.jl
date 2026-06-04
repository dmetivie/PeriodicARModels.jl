using Dates, Distributions, LinearAlgebra, Random


#### Simulating scenarios ####


function SimulateScenario!(L, p, n2t, Φ, σ::AbstractVector, rng)
    for (i, t) in enumerate(n2t[p+1:end])
        L[i+p] = sum(Φ[t, j] * L[i+p-j] for j in 1:p) + rand(rng, Normal(0, σ[t]))
    end
    return L
end

"""
    SimulateScenarios(x0::AbstractVector, Date_vec::AbstractVector, Φ, σ, nspart=0 ; n::Integer=1)

Simulate n temperature scenarios during the Date_vec timeline following an AR model with parameters Φ and σ and non stationnary part (trend + periodicity) nspart.
"""
function SimulateScenario(x0::AbstractArray, n2t::AbstractVector, Φ, σ::AbstractVector{T}, nspart, σ_nspart; rng=Random.default_rng(), correction="null", return_res=false) where T<:AbstractFloat
    p = length(x0)
    L0 = zeros(eltype(x0), length(n2t))
    L0[1:p] = x0
    L = SimulateScenario!(copy(L0), p, n2t, Φ, σ, rng)
    return return_res ? (L .* σ_nspart .+ nspart, L) : L .* σ_nspart .+ nspart
end
function SimulateScenario(x0::AbstractArray, n2t::AbstractVector, Φ, σ::AbstractVector{T}; rng=Random.default_rng(), correction="null", return_res=false) where T<:AbstractFloat
    p = length(x0)
    L0 = zeros(eltype(x0), length(n2t))
    L0[1:p] = x0
    L = SimulateScenario!(copy(L0), p, n2t, Φ, σ, rng)
    return L
end



#### Simulating paired scenarios ####


function SimulatePairedScenario!(L, p, n2t, Φ, Σ, d, rng)
    for (i, t) in enumerate(n2t[p+1:end])
        L[i+p, :] = sum(Φ[t][j] * L[i+p-j, :] for j in 1:p) .+ Σ[t] * randn(rng, d)
    end
    return L
end


#For the resample correction
function SimulatePairedScenarioCor!(L, p, n2t, Φ, Σ, d, σ_nspart_, nspart_, rng; max_resample=1000)
    for (i, t) in enumerate(n2t[p+1:end])
        L[i+p, :] = sum(Φ[t][j] * L[i+p-j, :] for j in 1:p) .+ Σ[t] * randn(rng, d)
        n_try = 1
        while @views L[i+p, 1] * σ_nspart_[i+p, 1] + nspart_[i+p, 1] > L[i+p, 2] * σ_nspart_[i+p, 2] + nspart_[i+p, 2] #While it is not good it tries again
            L[i+p, :] = sum(Φ[t][j] * L[i+p-j, :] for j in 1:p) .+ Σ[t] * randn(rng, d)
            n_try += 1
            n_try >= max_resample && (@warn "Max resampling ($max_resample) reached at step $(i+p); ordering constraint may not hold."; break)
        end
    end
    return L
end
function SimulatePairedScenarioCor!(L, p, n2t, Φ, Σ, d, rng; max_resample=1000)
    for (i, t) in enumerate(n2t[p+1:end])
        L[i+p, :] = sum(Φ[t][j] * L[i+p-j, :] for j in 1:p) .+ Σ[t] * randn(rng, d)
        n_try = 1
        while @views L[i+p, 1] > L[i+p, 2] #While it is not good it tries again
            L[i+p, :] = sum(Φ[t][j] * L[i+p-j, :] for j in 1:p) .+ Σ[t] * randn(rng, d)
            n_try += 1
            n_try >= max_resample && (@warn "Max resampling ($max_resample) reached at step $(i+p); ordering constraint may not hold."; break)
        end
    end
    return L
end


function SimulateScenario(x0::AbstractArray, n2t::AbstractVector, Φ, Σ::AbstractVector{Matrix{T}}, nspart, σ_nspart; rng=Random.default_rng(), correction="null", return_res=false) where T<:AbstractFloat
    p, d = collectpdx0(x0)
    M0 = copy(x0)
    p == 1 ? M0 = reshape(M0, (1, d)) : nothing
    L0 = zeros(eltype(x0), (length(n2t), d))
    L0[1:p, :] = M0

    L = correction == "resample" ? SimulatePairedScenarioCor!(copy(L0), p, n2t, Φ, Σ, d, σ_nspart, nspart, rng) : SimulatePairedScenario!(copy(L0), p, n2t, Φ, Σ, d, rng)

    return return_res ? ((L .* σ_nspart .+ nspart), L) : L .* σ_nspart .+ nspart
end
function SimulateScenario(x0::AbstractArray, n2t::AbstractVector, Φ, Σ::AbstractVector{Matrix{T}}; rng=Random.default_rng(), correction="null", return_res=false) where T<:AbstractFloat
    p, d = collectpdx0(x0)
    M0 = copy(x0)
    p == 1 ? M0 = reshape(M0, (1, d)) : nothing
    L0 = zeros(eltype(x0), (length(n2t), d))
    L0[1:p, :] = M0

    L = correction == "resample" ? SimulatePairedScenarioCor!(copy(L0), p, n2t, Φ, Σ, d, rng) : SimulatePairedScenario!(copy(L0), p, n2t, Φ, Σ, d, rng)

    return L
end

# SimulateScenarios(x0, Date_vec::AbstractVector, Φ, σ, nspart=0, rng=Random.default_rng(); n::Integer=1, index_nspart=nothing) = [SimulateScenario(x0, Date_vec, Φ, σ, nspart, rng, index_nspart=index_nspart) for _ in 1:n]




"""
    concat2by2(L)

(DEPRECATED)
Return the vector of the index-wise concatanations of the nested sub-sub-vectors in L.
The sub-vectors (not necessarily sub-sub-vectors) in L must have the same size. 
For example : 
 ```julia
    x = [[a₁, a₂],[a₃]]
    y = [[b₁], [b₂, b₃]]
    z = [[c₁], [c₂]]
    concat2by2([x,y,z]) = [[a₁, a₂, b₁, c₁], [a₃, b₂, b₃, c₂]]
```
"""
concat2by2(L1::AbstractVector, L2::AbstractVector) = [[u; v] for (u, v) in zip(L1, L2)]
concat2by2(L::AbstractVector) = reduce(concat2by2, L)

#Except this one
concat2by2mat(L1::AbstractVector, L2::AbstractVector) = [[u v] for (u, v) in zip(L1, L2)]
concat2by2mat(L::AbstractVector) = reduce(concat2by2mat, L)