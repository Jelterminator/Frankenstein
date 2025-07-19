module Analysis

using ..Core
using ..Utilities.Jacobians
using SciMLBase
using SparseArrays, Symbolics, ForwardDiff, LinearAlgebra

# Include submodule files
include("sparsity_analysis.jl")
include("stiffness_analysis.jl")
include("timescale_analysis.jl")
include("coupling_analysis.jl")
include("condition_analysis.jl")

# Export functions
export analyze_system_structure, detect_sparsity_patterns
export initial_stiffness_estimate, update_stiffness!
export compute_timescales, update_timescales!
export compute_coupling_strength, update_coupling_strength!
export compute_condition_number, update_condition_number!
export needs_analysis_update!

"""
    analyze_system_structure(prob::SciMLBase.ODEProblem) -> SystemAnalysis

Performs an initial comprehensive analysis of the ODE problem to populate a `SystemAnalysis` struct
with properties such as stiffness ratio, sparsity pattern, timescales, coupling strength, and
condition number. This function is computationally expensive and intended for use before integration
to guide solver configuration.

### Arguments
- `prob`: The ODE problem (`SciMLBase.ODEProblem`) to analyze, containing the function `f`, initial
  condition `u0`, parameters `p`, and time span `tspan`.

### Returns
- A `SystemAnalysis` struct containing:
  - `stiffness_ratio`: Estimated stiffness score.
  - `sparsity_pattern`: Jacobian sparsity pattern or `nothing`.
  - `timescales`: Vector of system timescales.
  - `coupling_strength`: Measure of variable coupling.
  - `condition_number`: Jacobian condition number.
  - `system_size`: Number of variables.
  - `is_sparse`: Boolean indicating if the system is sparse.
  - `jacobian`: Cached Jacobian matrix (sparse if applicable).
  - `stable_count`, `last_update_step`, `current_step`, `last_norm_du`, `history`: State tracking fields.
"""

function analyze_system_structure(prob)
    u0 = prob.u0
    t0 = prob.tspan[1]
    p = prob.p
    f = prob.f

    # Compute system size and sparsity
    system_size = length(u0)
    sparsity = detect_sparsity_patterns(prob)
    is_sparse = sparsity !== nothing && nnz(sparsity) / (system_size^2) < 0.1

    # Compute Jacobian with fallback
    try
        J = compute_jacobian(f, u0, p, t0)
    catch e
        @warn "Jacobian computation failed: $e. Using finite differences."
        J = finite_difference_jac(u -> SciMLBase.isinplace(f) ? (du = similar(u); f(du, u, p, t0); du) : f(u, p, t0), u0, p)
    end
    J = is_sparse ? sparse(J) : J  # Convert to sparse if applicable

    # Perform analyses
    stiffness = Stiffness.initial_stiffness_estimate(f, u0, p, J0=J)
    timescales = compute_timescales(prob, u0, t0, J=J)
    coupling = compute_coupling_strength(prob, u0, t0, J=J)
    condition = compute_condition_number(prob, u0, t0, J=J)

    return SystemAnalysis{Float64}(stiffness, sparsity, timescales, coupling, condition, system_size, is_sparse, J, 0)
end

"""
    needs_analysis_update!(sa::SystemAnalysis, step_info::Core.StepInfo; 
        step_change_thresh=2.0, reject_thresh=0.2, norm_change_thresh=0.5, 
        stable_steps=10, min_update_interval=5) -> NamedTuple

Determines which analysis variables need updating based on dynamic indicators.
Returns a NamedTuple with boolean flags indicating which updates are required.

### Arguments
- `sa`: The `SystemAnalysis` struct containing current system properties and state.
- `step_info`: A `StepInfo` struct with integrator data (`u`, `du`, `dt`, `dt_prev`, `rejects`).
- `step_change_thresh`: Threshold for step size change ratio to trigger updates.
- `reject_thresh`: Threshold for rejection rate to trigger updates.
- `norm_change_thresh`: Threshold for solution norm change to trigger updates.
- `stable_steps`: Number of stable steps before reducing update frequency.
- `min_update_interval`: Minimum number of steps between updates to avoid over-triggering.

### Returns
- A NamedTuple with fields:
  - `stiffness`: Bool, whether to update `stiffness_ratio`.
  - `timescales`: Bool, whether to update `timescales`.
  - `coupling`: Bool, whether to update `coupling_strength`.
  - `condition`: Bool, whether to update `condition_number`.
"""
function needs_analysis_update!(sa::SystemAnalysis, step_info::Core.StepInfo; 
        step_change_thresh=2.0, reject_thresh=0.2, norm_change_thresh=0.5, 
        stable_steps=10, min_update_interval=5)
    # Initialize update flags
    update_stiffness = false
    update_timescales = false
    update_coupling = false
    update_condition = false

    # Get or initialize last update step
    last_update_step = get!(sa, :last_update_step, 0)
    current_step = get!(sa, :current_step, 0) + 1
    sa.current_step = current_step

    # Skip updates if within minimum interval, unless forced by significant changes
    if current_step - last_update_step < min_update_interval && sa.stable_count < stable_steps
        return (stiffness=false, timescales=false, coupling=false, condition=false)
    end

    # Compute lightweight indicators
    dt_ratio = step_info.dt_prev > 0 ? step_info.dt_prev / step_info.dt : 1.0
    reject_rate = step_info.rejects / max(1, get(sa, :history, 10))
    norm_du = norm(step_info.du)
    norm_u = norm(step_info.u)
    norm_change = get!(sa, :last_norm_du, norm_du) > 0 ? abs(norm_du / sa.last_norm_du - 1) : 0.0

    # Update last norm for next call
    sa.last_norm_du = norm_du

    # Check conditions for updates
    significant_change = dt_ratio > step_change_thresh || reject_rate > reject_thresh || norm_change > norm_change_thresh

    if significant_change || sa.stable_count >= stable_steps
        # Stiffness: Trigger on step size changes or high rejections
        update_stiffness = dt_ratio > step_change_thresh || reject_rate > reject_thresh

        # Timescales: Trigger on solution norm changes or significant step changes
        update_timescales = norm_change > norm_change_thresh || dt_ratio > step_change_thresh

        # Coupling: Trigger on solution norm changes or high rejections
        update_coupling = norm_change > norm_change_thresh || reject_rate > reject_thresh

        # Condition: Trigger on step size changes or high rejections
        update_condition = dt_ratio > step_change_thresh || reject_rate > reject_thresh

        # Reset stable count if significant change detected
        if significant_change
            sa.stable_count = 0
        end

        # Update last update step if any update is triggered
        if update_stiffness || update_timescales || update_coupling || update_condition
            sa.last_update_step = current_step
        end
    end

    return (stiffness=update_stiffness, timescales=update_timescales, 
            coupling=update_coupling, condition=update_condition)
end

end # module Analysis