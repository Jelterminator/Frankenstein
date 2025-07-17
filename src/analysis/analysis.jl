# Analysis.jl

module Analysis

using ..Core
using SciMLBase
using SparseArrays, Symbolics, ForwardDiff, LinearAlgebra

# Include submodule files
include("sparsity_analysis.jl")
include("stiffness_analysis.jl")
include("timescale_analysis.jl")
include("coupling_analysis.jl")
include("condition_analysis.jl")

# Export functions
export analyze_system_structure
export detect_sparsity_patterns
export initial_stiffness_estimate, update_stiffness!
export compute_timescales, compute_coupling_strength
export compute_condition_number

"""
    analyze_system_structure(prob) -> SystemAnalysis

Analyze the ODE problem to populate the SystemAnalysis struct with stiffness, sparsity,
timescales, coupling strength, condition number, and Jacobian matrix.
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

    # Compute Jacobian once for all analyses
    J = ConditionAnalysis.compute_jacobian(f, u0, p, t0)

    # Perform analyses using cached Jacobian
    stiffness = Stiffness.initial_stiffness_estimate(prob, J=J)
    timescales = compute_timescales(prob, u0, t0, J=J)
    coupling = compute_coupling_strength(prob, u0, t0, J=J)
    condition = compute_condition_number(prob, u0, t0, J=J)

    return SystemAnalysis{Float64}(stiffness, sparsity, timescales, coupling, condition, system_size, is_sparse, J)
end

end # module Analysis