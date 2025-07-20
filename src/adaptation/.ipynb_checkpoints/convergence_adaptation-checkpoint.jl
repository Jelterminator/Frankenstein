## adaptation/convergence_adaptation.jl
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

function adapt!(strategy::ConvergenceAdaptation, analysis::SystemAnalysis, rec::AlgorithmRecommendation)
    # 1. Evaluate recent step errors
    # 2. Increase/decrease order or switch method family
    return rec
end

export ConvergenceAdaptation, adapt!
end