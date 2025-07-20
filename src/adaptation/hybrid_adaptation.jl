## adaptation/hybrid_adaptation.jl
module HybridAdaptation

using ..Core: SystemAnalysis
using ..Core: AbstractAdaptationStrategy
using ..Core: AlgorithmRecommendation

"""
    HybridAdaptation

Combine multiple adaptation mechanisms: solver-switching, backend-switching, preconditioner retuning.
"""
struct HybridAdaptation <: AbstractAdaptationStrategy
    # composite strategy holding sub-strategies
    performance::PerformanceAdaptation
    stability::StabilityAdaptation
    convergence::ConvergenceAdaptation
    memory::MemoryAdaptation
    parallel::ParallelAdaptation
end

function adapt!(strategy::HybridAdaptation, analysis::SystemAnalysis, rec::AlgorithmRecommendation)
    # apply each sub-strategy in sequence
    rec = adapt!(strategy.performance, analysis, rec)
    rec = adapt!(strategy.stability, analysis, rec)
    rec = adapt!(strategy.convergence, analysis, rec)
    rec = adapt!(strategy.memory, analysis, rec)
    rec = adapt!(strategy.parallel, analysis, rec)
    return rec
end

export HybridAdaptation, adapt!
end