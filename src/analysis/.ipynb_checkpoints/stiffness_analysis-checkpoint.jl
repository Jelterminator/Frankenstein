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
    update_stiffness!(sa::SystemAnalysis, step_info; 
        dyn_weights=(1.0,1.0,1.0), dyn_max=log10(101), thresh=100, history=10,
        decay_rate=0.9, stable_steps=5
    )

Update the SystemAnalysis.sa.stiffness_ratio based on dynamic step behavior,
trigger deeper analysis as needed, and apply decay when system is stable:
 1. Step shrink factor
 2. Reject ratio
 3. Dynamic derivative ratio
 4. Decay stiffness when indicators show non-stiff behavior
"""
function update_stiffness!(sa::SystemAnalysis, step_info; 
        dyn_weights=(1.0,1.0,1.0), dyn_max=log10(101), thresh=100, history=10,
        decay_rate=0.9, stable_steps=5)
    # Dynamic indicators
    shrink = max(0.0, step_info.dt_prev > 0 ? log10(step_info.dt_prev/step_info.dt) : 0.0)
    rejects = step_info.rejects / history
    du = step_info.du
    Ddyn = maximum(abs.(du)) / max(1e-12, minimum(abs.(du[abs.(du) .> 0])))

    # Raw dynamic score
    raw = dyn_weights[1]*shrink + dyn_weights[2]*log10(1+rejects) + dyn_weights[3]*log10(1+Ddyn)
    dyn_score = clamp(raw / dyn_max * 10000, 0, 10000)

    # Decide if system is showing non-stiff behavior
    if dyn_score < thresh/2
        sa.stiffness_ratio *= decay_rate
        sa.stable_count = get(sa, :stable_count, 0) + 1
    else
        sa.stable_count = 0
        # Combine with previous
        sa.stiffness_ratio = clamp((sa.stiffness_ratio + dyn_score)/2, 0, 10000)
    end

    sa.is_stiff = sa.stiffness_ratio >= thresh

    return nothing
end

end # module Analysis.Stiffness
