module Analysis.Stiffness

using LinearAlgebra, ForwardDiff
using ..Core: SystemAnalysis

export initial_stiffness_estimate, update_stiffness!

"""
    gershgorin_spectral_bound(J::AbstractMatrix) -> Float64

Compute an upper bound on the spectral radius of J using Gershgorin circles.
"""
function gershgorin_spectral_bound(J::AbstractMatrix)
    diags = diag(J)
    row_sums = sum(abs.(J), dims=2) .- abs.(diags)
    return maximum(abs.(diags) .+ row_sums)
end

"""
    initial_stiffness_estimate(f, u0, p; 
        J_method = :fd, ngersh = 3, weights = (1.0,1.0,1.0,0.0), raw_max = 12.0
    ) -> Float64

Compute a pre-integration stiffness score in [0,10000] based on:
 1. Gershgorin spectral bound
 2. Derivative ratio at u0
 3. Jacobian norm at u0
 4. Optional user hint (via weights[4])
"""
# stiffness_analysis.jl (partial update)

function initial_stiffness_estimate(f, u0, p; J_method=:fd, ngersh=3, weights=(1.0,1.0,1.0,0.0), raw_max=12.0, J0=nothing)
    # Use provided J0 for initial Jacobian if available
    if J0 === nothing
        J0 = J_method == :fd ? finite_difference_jac(f, u0, p) :
              ForwardDiff.jacobian(u -> f(u, p), u0)
    end
    # Compute Gershgorin bounds at perturbed points (still needs multiple Jacobians)
    Gs = Float64[]
    for i in 1:ngersh
        u = u0 .+ 1e-6 * randn(length(u0))
        J = J_method == :fd ? finite_difference_jac(f, u, p) :
            ForwardDiff.jacobian(u -> f(u, p), u)
        push!(Gs, gershgorin_spectral_bound(J))
    end
    G = maximum(Gs)

    # Use J0 for Jacobian norm
    Jn = opnorm(J0, Inf)

    # Rest of the computation...
    du = f(u0, p)
    D = maximum(abs.(du)) / max(1e-12, minimum(abs.(du[abs.(du) .> 0])))
    raw = weights[1]*log10(1+G) + weights[2]*log10(1+D) + weights[3]*log10(1+Jn) + weights[4]*0
    score = clamp(raw / raw_max * 10000, 0, 10000)
    return score
end

"""
    update_stiffness!(sa::SystemAnalysis, step_info::Core.StepInfo; 
        dyn_weights=(0.5,0.5,0.5,0.5), dyn_max=log10(101), thresh=100, history=10,
        decay_rate=0.9, stable_steps=5, recompute_jacobian_interval=20)

Updates the `stiffness_ratio` in `SystemAnalysis` based on dynamic indicators (step size changes,
rejection rates, derivative ratios, and Gershgorin spectral bound). Designed to be computationally
cheap by reusing the cached Jacobian when possible and limiting Jacobian recomputations.

### Arguments
- `sa`: The `SystemAnalysis` struct containing current system properties.
- `step_info`: A `StepInfo` struct with integrator data (`u`, `du`, `dt`, `dt_prev`, `rejects`).
- `dyn_weights`: Weights for combining step shrink, rejections, derivative ratio, and Gershgorin bound.
- `dyn_max`: Normalization factor for dynamic score.
- `thresh`: Stiffness threshold for setting `is_stiff`.
- `history`: Number of steps to consider for rejection rate.
- `decay_rate`: Factor to reduce stiffness when non-stiff behavior is detected.
- `stable_steps`: Number of stable steps before applying decay.
- `recompute_jacobian_interval`: Minimum steps between Jacobian recomputations.

### Effects
- Updates `sa.stiffness_ratio`, `sa.is_stiff`, and `sa.stable_count` in-place.
- Recomputes `sa.jacobian` if necessary and allowed.
"""
function update_stiffness!(sa::SystemAnalysis, step_info::Core.StepInfo; 
        dyn_weights=(0.5,0.5,0.5,0.5), dyn_max=log10(101), thresh=100, history=10,
        decay_rate=0.9, stable_steps=5, recompute_jacobian_interval=20)
    # Dynamic indicators
    shrink = max(0.0, step_info.dt_prev > 0 ? log10(step_info.dt_prev/step_info.dt) : 0.0)
    rejects = step_info.rejects / history
    du = step_info.du
    Ddyn = maximum(abs.(du)) / max(1e-12, minimum(abs.(du[abs.(du) .> 0])))

    # Gershgorin spectral bound
    G = 0.0
    if sa.jacobian !== nothing
        G = gershgorin_spectral_bound(sa.jacobian)
    end

    # Recompute Jacobian if significant change detected and interval elapsed
    if sa.current_step - get(sa, :last_jacobian_update, 0) >= recompute_jacobian_interval
        update_flags = Analysis.needs_analysis_update!(sa, step_info)
        if update_flags.stiffness && sa.jacobian !== nothing
            try
                sa.jacobian = Utilities.Jacobians.compute_jacobian(
                    step_info.prob.f, step_info.u, step_info.prob.p, step_info.t)
                sa.jacobian = sa.is_sparse ? sparse(sa.jacobian) : sa.jacobian
                sa.last_jacobian_update = sa.current_step
                G = gershgorin_spectral_bound(sa.jacobian)
            catch e
                @warn "Jacobian recomputation failed: $e. Using cached Jacobian."
            end
        end
    end

    # Combine indicators
    raw = (dyn_weights[1]*shrink + dyn_weights[2]*log10(1+rejects) + 
           dyn_weights[3]*log10(1+Ddyn) + dyn_weights[4]*log10(1+G))
    dyn_score = clamp(raw / dyn_max * 10000, 0, 10000)

    # Update stiffness with decay for non-stiff behavior
    if dyn_score < thresh/2
        sa.stiffness_ratio *= decay_rate
        sa.stable_count += 1
    else
        sa.stable_count = 0
        sa.stiffness_ratio = clamp((sa.stiffness_ratio + dyn_score)/2, 0, 10000)
    end

    sa.is_stiff = sa.stiffness_ratio >= thresh
    return nothing
end

end # module Analysis.Stiffness
