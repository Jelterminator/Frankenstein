# condition_analysis.jl

function compute_condition_number(prob, u=prob.u0, t=prob.tspan[1]; J=nothing)
    if J === nothing
        J = Utilities.Jacobians.compute_jacobian(prob.f, u, prob.p, t)
    end
    return cond(J)
end

function update_condition_number!(sa::SystemAnalysis, step_info::Core.StepInfo)
    # Approximate condition number change based on stiffness and du
    du = step_info.du
    norm_du = norm(du)
    norm_u = norm(step_info.u)
    est_cond = sa.condition_number * (norm_du / max(norm_u, 1e-10))  # Heuristic scaling
    sa.condition_number = clamp(0.5 * (sa.condition_number + est_cond), 1.0, 1e10)
    return nothing
end