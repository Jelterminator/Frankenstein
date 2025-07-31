module StabilityAdaptation

using ..Core: SystemAnalysis
using ..Core: AbstractAdaptationStrategy
using ..Core: AlgorithmRecommendation

"""
    StabilityAdaptation

Switch solver methods or parameters to maintain numerical stability under dynamic stiffness.
"""
struct StabilityAdaptation <: AbstractAdaptationStrategy
    stiffness_threshold::Float64
end

"""
    adapt!(strategy::StabilityAdaptation, analysis::SystemAnalysis, rec::AlgorithmRecommendation)

Monitor system stiffness and update the algorithm recommendation if stability is at risk.
"""
function adapt!(strategy::StabilityAdaptation, analysis::SystemAnalysis, rec::AlgorithmRecommendation)
    # 1. Estimate current stiffness metric from analysis
    local_stiffness = getfield(analysis, :stiffness_estimate)

    # 2. If stiffness exceeds threshold, switch to a more stable implicit solver
    if local_stiffness > strategy.stiffness_threshold
        # Log or mark the adaptation
        println("[StabilityAdaptation] High stiffness detected: ", local_stiffness,
                " > ", strategy.stiffness_threshold)

        # Create a new recommendation: 
        new_algo =select_best_algorithm(analysis)  
        # Preserve existing tolerance and other parameters if present
        params = hasproperty(rec, :params) ? rec.params : Dict{Symbol,Any}()
        params[:method] = new_algo

        return AlgorithmRecommendation(new_algo; params...)
    end

    # Otherwise, leave recommendation unchanged
    return rec
end

end # module StabilityAdaptation