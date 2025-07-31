module ConvergenceAdaptation

using ..Core: SystemAnalysis
using ..Core: AbstractAdaptationStrategy
using ..Core: AlgorithmRecommendation

"""
    ConvergenceAdaptation

Adjust step-size and solver order to ensure convergence targets are met.
"""
struct ConvergenceAdaptation <: AbstractAdaptationStrategy
    error_tolerance::Float64
end

"""
    adapt!(strategy::ConvergenceAdaptation, analysis::SystemAnalysis, rec::AlgorithmRecommendation)

Monitor recent convergence errors and adjust solver parameters accordingly.
"""
function adapt!(strategy::ConvergenceAdaptation, analysis::SystemAnalysis, rec::AlgorithmRecommendation)
    # 1. Extract current error estimate from analysis
    local_error = getfield(analysis, :error_estimate)

    # Prepare parameters dict
    params = hasproperty(rec, :params) ? copy(rec.params) : Dict{Symbol,Any}()

    # 2. If error too large, increase solver order or tighten tolerance
    if local_error > strategy.error_tolerance
        println("[ConvergenceAdaptation] Error ", local_error,
                " exceeds tolerance ", strategy.error_tolerance)
        # Increase solver order if supported
        current_order = get(params, :order, 1)
        params[:order] = current_order + 1
    # 3. If error is much smaller, consider reducing solver order
    elseif local_error < strategy.error_tolerance / 10
        println("[ConvergenceAdaptation] Error ", local_error,
                " well below tolerance ", strategy.error_tolerance)
        current_order = get(params, :order, 1)
        params[:order] = max(1, current_order - 1)
    end

    # Return updated recommendation
    return AlgorithmRecommendation(rec.method; params...)
end

export ConvergenceAdaptation, adapt!
end # module ConvergenceAdaptation
