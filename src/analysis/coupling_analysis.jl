# coupling_analysis.jl

using LinearAlgebra

function compute_coupling_strength(prob, u=prob.u0, t=prob.tspan[1]; J=nothing)
    if J === nothing
        J = ConditionAnalysis.compute_jacobian(prob.f, u, prob.p, t)
    end
    n = size(J, 1)
    sum_diag = sum(abs(J[i,i]) for i in 1:n)
    sum_off_diag = sum(abs(J[i,j]) for i in 1:n for j in 1:n if i != j)
    coupling_strength = sum_off_diag / max(sum_diag, 1e-10)
    return coupling_strength
end

function update_coupling_strength!(sa::SystemAnalysis, step_info::Core.StepInfo)
    du = step_info.du
    # Estimate coupling from correlation of component changes
    n = length(du)
    sum_off_diag = sum(abs(du[i] - du[j]) for i in 1:n for j in i+1:n) / (n*(n-1)/2)
    sum_diag = sum(abs.(du)) / n
    new_coupling = sum_off_diag / max(sum_diag, 1e-10)
    sa.coupling_strength = 0.5 * (sa.coupling_strength + new_coupling)  # Smooth update
    return nothing
end