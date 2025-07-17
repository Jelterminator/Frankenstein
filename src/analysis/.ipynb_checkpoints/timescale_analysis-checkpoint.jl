# timescale_analysis.jl

using LinearAlgebra

function compute_timescales(prob, u=prob.u0, t=prob.tspan[1]; J=nothing)
    if J === nothing
        J = ConditionAnalysis.compute_jacobian(prob.f, u, prob.p, t)
    end
    evs = eigvals(J)
    reals = real.(evs)
    timescales = 1.0 ./ max.(abs.(reals), 1e-10)
    return timescales
end